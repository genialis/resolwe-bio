"""Call peaks with MACS 2.0."""

import os
import re
from pathlib import Path

import pandas as pd
from pandas.errors import EmptyDataError
from plumbum import RETCODE, TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    GroupField,
    IntegerField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)

BED_COLUMNS = ["chromosome", "start", "end", "name", "score", "strand"]

BEDPE_COLUMNS = [
    "chromosome",
    "start",
    "end",
    "chromosome_r2",
    "start_r2",
    "end_r2",
    "name",
    "score",
    "strand",
    "strand_r2",
]

SPECIES_GSIZES = {
    "Homo sapiens": "2.7e9",
    "Mus musculus": "1.87e9",
    "Dictyostelium discoideum": "3.4e7",
    "Drosophila melanogaster": "1.2e8",
    "Caenorhabditis elegans": "9e7",
    "Saccharomyces cerevisiae": "12.1e6",
    "Rattus norvegicus": "2e9",
}

SPP_HEADER = [
    "Reads",
    "Est. Fragment Len.",
    "Corr. Est. Fragment Len.",
    "Phantom Peak",
    "Corr. Phantom Peak",
    "Argmin. Corr.",
    "Min. Corr.",
    "NSC",
    "RSC",
    "Quality Tag",
]


def test_paired_end(bam_file, error):
    """Test if a given bam file is paired end."""
    args = ["-c", "-f", "1", bam_file]
    return_code, n_paired, stderr = Cmd["samtools"]["view"][args] & TEE(retcode=None)
    if return_code:
        print(stderr)
        error("Failed to count the number of paired end reads.")

    try:
        n_paired = int(n_paired)
    except ValueError:
        error("Could not determine if the reads are single or paired end.")

    return n_paired > 0


def merge_dict(d1, d2):
    """Merge two dictionaries."""
    return {**d1, **d2}


def filter_bam(bam_file, out_bam, min_quality, is_paired, name, error):
    """Filter bam file."""
    # Remove unmapped reads, not primary alignments, duplicates.
    filter_params = [
        "-F",
        1804,
        "-b",
        "-o",
        out_bam,
    ]

    # Remove reads below minimum mapping quality.
    if min_quality > 0:
        filter_params.extend(["-q", min_quality])

    # Remove reads not mapped in proper pair.
    if is_paired:
        filter_params.extend(["-f", 2])

    filter_params.append(bam_file)

    return_code, stdout, stderr = Cmd["samtools"]["view"][filter_params] & TEE(
        retcode=None
    )

    if return_code:
        print(stdout, stderr)
        error(f"Samtools filtering failed for {name}.")


def name_sort_bam(in_bam, out_bam, error):
    """Name sort bam file."""
    sort_params = ["-n", "-o", out_bam, in_bam]
    return_code, stdout, stderr = Cmd["samtools"]["sort"][sort_params] & TEE(
        retcode=None
    )

    if return_code:
        print(stdout, stderr)
        error(f"Samtools sorting failed for {in_bam}.")


def is_empty(path):
    """Check if a certain file is empty."""
    return Path(path).is_file() and Path(path).stat().st_size > 0


def parse_markdup(markdup_stats):
    """Parse mark duplicates report."""
    with open(markdup_stats, "r") as handle:
        lines = handle.readlines()
        if "## HISTOGRAM\tjava.lang.Double\n" in lines:
            histogram_start = lines.index("## HISTOGRAM\tjava.lang.Double\n")
            skip_lines = list(range(histogram_start, len(lines)))
        else:
            skip_lines = None

    return (
        pd.read_csv(markdup_stats, sep="\t", comment="#", skiprows=skip_lines)
        .squeeze()
        .to_dict()
    )


def read_bed(bed_file, col_names):
    """Read a BED file and add a header."""
    return pd.read_csv(
        bed_file,
        sep="\t",
        header=None,
        names=col_names,
        float_precision="round_trip",
    )


def drop_mt(bed):
    """Drop mitochondrial genome from the bed file."""
    bed = bed.drop(bed[bed.chromosome.isin(["chrM", "MT"])].index)
    if "chromosome_r2" in bed.columns:
        bed = bed.drop(bed[bed.chromosome.isin(["chrM", "MT"])].index)
    return bed


def parse_flagstat(report, warning):
    """Parse flagstat report and extract read statistics."""

    # Typical lines we are interested in would be:
    # 400 + 0 in total (QC-passed reads + QC-failed reads)
    # 396 + 0 mapped (99.00% : N/A)
    flagstat_regexes = {
        "TOTAL_READS": r"(\d+) \+ \d+ in total \(QC-passed reads \+ QC-failed reads\)",
        "MAPPED_READS": r"(\d+) \+ \d+ mapped \(.+:.+\)",
        "MAPPED_PERCENTAGE": r"\d+ \+ \d+ mapped \((.+):.+\)",
    }

    parsed_data = {}
    for k, r in flagstat_regexes.items():
        r_search = re.search(r, report, re.MULTILINE)
        try:
            parsed_data[k] = r_search.group(1).strip()
        except (AttributeError, IndexError):
            # No matches found.
            parsed_data[k] = "NA"
            warning(f"Failed to determine {k} based on flagstat report.")
    return parsed_data


def parse_bowtie2_report(report_path, n_mapped, warning):
    """Parse Bowtie2 alignment report."""
    with open(report_path, "r") as handle:
        report = handle.read()

    # The report is structured as:
    # 1000000 reads; of these:
    #   1000000 (100.00%) were paired; of these:
    #        980354 (98.04%) aligned concordantly 0 times
    #   ...
    # 5.53% overall alignment rate
    stats_regexes = {
        "TOTAL_READS": r"(\d+) reads; of these:",
        "MAPPED_PERCENTAGE": r"(\d+\.\d+\%) overall alignment rate",
    }

    parsed_data = {
        "TOTAL_READS": "NA",
        "MAPPED_READS": n_mapped,
        "MAPPED_PERCENTAGE": "NA",
    }
    for k, r in stats_regexes.items():
        r_search = re.search(r, report, re.MULTILINE)
        try:
            parsed_data[k] = r_search.group(1).strip()
        except (AttributeError, IndexError):
            # No matches found.
            warning(f"Failed to determine {k} based on bowtie2 report.")

    if parsed_data["TOTAL_READS"].isdigit():
        parsed_data["TOTAL_READS"] = int(parsed_data["TOTAL_READS"]) * 2
    return parsed_data


def convert_bam(bam_file, name, out_bed, is_paired, error, first_mate=False):
    """Convert bam file to bed file."""

    args = ["-i", bam_file]

    if is_paired:
        args.append("-bedpe")
        if first_mate:
            args.append("-mate1")

    return_code = (Cmd["bedtools"]["bamtobed"][args] > out_bed) & RETCODE
    if return_code:
        error(f"Failed to convert bam file to bed file for {name}.")


def get_pbc_metrics(bam_file, name, out_bed, is_paired, error):
    """Get PCR bottlenecking coefficient metrics.

    Metrics are the following:
    Total Reads: number of reads
    Distinct Reads: number of all genomic locations where reads mapped
    One Read: number of genomic locations where only one read maps uniquely
    Two Reads: number of genomic locations where 2 reads map uniquely
    NRF = Non-Redundant Fraction (Distinct Reads / Total Reads)
    PBC1 = PCR Bottlenecking Coefficient 1 (One Read / Distinct Reads)
    PBC2 = PCR Bottlenecking Coefficient 2 (One Read / Two Reads)

    Matching bash scripts available on the link below.
    https://github.com/ENCODE-DCC/chip-seq-pipeline/blob/master/dnanexus/filter_qc/src/filter_qc.py
    """

    convert_bam(
        bam_file=bam_file, name=name, out_bed=out_bed, is_paired=is_paired, error=error
    )

    if is_paired:
        col_names = BEDPE_COLUMNS
        column_subset = [
            "chromosome",
            "start",
            "chromosome_r2",
            "end_r2",
            "strand",
            "strand_r2",
        ]
    else:
        col_names = BED_COLUMNS
        column_subset = ["chromosome", "start", "end", "strand"]

    bed = read_bed(out_bed, col_names=col_names)

    bed = bed.loc[:, column_subset]
    bed = drop_mt(bed)
    # Add a column with the number of duplicate entires.
    bed = bed.groupby(column_subset).size().reset_index(name="count")

    pbc = {
        "Distinct Reads": len(bed.index),
        "One Read": sum(bed["count"] == 1),
        "Two Reads": sum(bed["count"] == 2),
        "NRF": "N/A",
        "PBC1": "N/A",
        "PBC2": "N/A",
        "Total Reads": bed["count"].sum(),
    }

    pbc["NRF"] = round(pbc["Distinct Reads"] / pbc["Total Reads"], 6)
    pbc["PBC1"] = round(pbc["One Read"] / pbc["Distinct Reads"], 6)
    if pbc["Two Reads"] > 0:
        pbc["PBC2"] = round(pbc["One Read"] / pbc["Two Reads"], 6)

    return pbc


def prepare_tagalign(bed_file, tagalign_file, is_paired, tn5_shift, compression=None):
    """Prepare a final tagAlign file.

    First transform deduplicated bam to tagAlign and then optionally do
    the Tn5 transposon shifting.
    """
    col_names = BED_COLUMNS
    if is_paired:
        tagalign = read_bed(bed_file, col_names=BEDPE_COLUMNS)

        r2_col_names = [
            "chromosome_r2",
            "start_r2",
            "end_r2",
            "name",
            "score",
            "strand_r2",
        ]

        name_mapper = dict(zip(r2_col_names, col_names))

        # Interleave R1 and R2 reads.
        tagalign = pd.concat(
            [
                tagalign[col_names].reset_index(drop=True),
                tagalign[r2_col_names]
                .rename(columns=name_mapper)
                .set_index(tagalign.index + 0.1),
            ],
            sort=False,
        )

        tagalign = tagalign.sort_index().reset_index(drop=True)

    else:
        tagalign = read_bed(
            bed_file,
            col_names=col_names,
        )

    # Remove sequence info in the name field and set scores to the max.
    tagalign["name"] = "N"
    tagalign["score"] = 1000

    if tn5_shift:
        tagalign.loc[(tagalign.strand == "+"), "start"] = tagalign.start + 4
        tagalign.loc[(tagalign.strand == "-"), "end"] = tagalign.end - 5

    tagalign.to_csv(
        tagalign_file,
        sep="\t",
        index=False,
        header=False,
        compression=compression,  # Pandas 0.23.4 compatibility
    )


def subsample_tagalign(tagalign_file, out_file, is_paired, n_sub, compression=None):
    """Subsample a tagAlign file and save it."""
    tagalign = read_bed(tagalign_file, col_names=BED_COLUMNS)

    if is_paired:
        # Subsample only R1 reads which are in every second row.
        tagalign = tagalign.iloc[::2, :]

    tagalign = drop_mt(tagalign)

    # Sample size has to be less or equal to the number of entries.
    n_sub = min(n_sub, len(tagalign.index))
    tagalign = tagalign.sample(n=n_sub, random_state=42)
    tagalign.to_csv(
        out_file, sep="\t", index=False, header=False, compression=compression
    )


def save_prepeak_qc(metrics, out_file):
    """Save a prepeak QC file."""
    prepeak_qc_fields = [
        "TOTAL_READS",
        "MAPPED_READS",
        "MAPPED_PERCENTAGE",
        "UNPAIRED_READS_EXAMINED",
        "READ_PAIRS_EXAMINED",
        "UNPAIRED_READ_DUPLICATES",
        "PERCENT_DUPLICATION",
        "NRF",
        "PBC1",
        "PBC2",
        "NSC",
        "RSC",
    ]

    pd.Series(metrics).loc[prepeak_qc_fields].to_frame().T.to_csv(
        out_file, sep="\t", index=False
    )


def correct_bed_file(in_bed, out_bed, error):
    """Correct bed file to be compatible with bedToBigBed tool.

    1.) restrict all scores to the maximal value of 1000,
    2.) in strand column replace '?' with '.'.
    """
    try:
        df = pd.read_csv(in_bed, delimiter="\t", header=None, dtype=str)
    except EmptyDataError:
        error(
            f"The input BED file {in_bed} is empty. Your analysis might "
            "have failed to identify regions of interest (peaks, junctions, etc.)."
        )

    df.iloc[:, 4] = pd.to_numeric(df.iloc[:, 4]).round().astype(int)
    df.iloc[:, 4] = df.iloc[:, 4].clip(upper=1000)

    # If strand column exist replace '?' with '.'.
    if len(df.columns) >= 6:
        df.iloc[:, 5] = df.iloc[:, 5].replace("?", ".")

    df.to_csv(out_bed, sep="\t", index=False, header=False)


def clip_bed(in_bed, out_bed, chromosome_sizes, file_type, error):
    """Clip bed file to remove off-chromosome places."""
    padded_bed = Path("padded.bed")
    slop_command = Cmd["bedtools"]["slop"][
        "-i", in_bed, "-g", chromosome_sizes, "-b", 0
    ] > str(padded_bed)

    return_code = slop_command & RETCODE

    if return_code:
        error(f"Failed to increase the size of features in {file_type} file.")

    return_code, _, stderr = Cmd["bedClip"][
        padded_bed, chromosome_sizes, out_bed
    ] & TEE(retcode=None)

    if return_code:
        print(stderr)
        error(
            "Preventing the extension of intervals beyond chromosome boundaries for "
            f"{file_type} file failed."
        )

    padded_bed.unlink()


def create_big_bed(in_file, out_bb, chromosome_sizes, file_type, error):
    """Create BigBed for IGV and UCSC genome browsers."""

    clipped_file = f"peaks_clip_{in_file}"
    clip_bed(
        in_bed=in_file,
        out_bed=clipped_file,
        chromosome_sizes=chromosome_sizes,
        file_type=f"*{file_type}",
        error=error,
    )

    if file_type == "narrowPeak":
        big_bed_params = [
            "as=/opt/kent/bedToBigBed/narrowPeak.as",
            "-type=bed6+4",
        ]
    else:
        big_bed_params = []

    big_bed_params.extend([clipped_file, chromosome_sizes, out_bb])

    return_code, stdout, stderr = Cmd["bedToBigBed"][big_bed_params] & TEE(retcode=None)

    if return_code:
        print(stdout, stderr)
        error(f"Creating BigBed from {file_type} file failed.")


def create_bigwig(in_bdg, out_bw, chromosome_sizes, file_type, error):
    """Create BigWig file from BedGraph file."""
    clipped_pileup_bdg = f"clipped_{in_bdg}"
    clip_bed(
        in_bed=in_bdg,
        out_bed=clipped_pileup_bdg,
        chromosome_sizes=chromosome_sizes,
        file_type=file_type,
        error=error,
    )

    old_locale = os.environ.get("LC_COLLATE", None)
    os.environ["LC_COLLATE"] = "C"
    sorted_bdg = f"sorted_{in_bdg}"
    sort_command = (
        Cmd["sort"]["-k", "1,1", "-k", "2,2n", clipped_pileup_bdg] > sorted_bdg
    )
    if old_locale:
        os.environ["LC_COLLATE"] = old_locale

    return_code = sort_command & RETCODE

    if return_code:
        error(f"Sorting the {in_bdg} file failed.")

    merged_bdg = f"sorted_merged_{in_bdg}"
    merge_command = (
        Cmd["bedtools"]["merge"]["-d", "-1", "-c", 4, "-o", "mean", "-i", sorted_bdg]
        > merged_bdg
    )

    return_code = merge_command & RETCODE

    if return_code:
        error(f"Interval merging for {in_bdg} file failed.")

    return_code, stdout, stderr = Cmd["bedGraphToBigWig"][
        merged_bdg, chromosome_sizes, out_bw
    ] & TEE(retcode=None)

    if return_code:
        print(stdout, stderr)
        error(f"Creating bigWig from {file_type} bedGraph failed")


def count_lines(path):
    """Count the number of lines in a file."""
    with open(path, "r") as f:
        return sum(1 for _ in f)


def count_overlap(file_a, file_b, error):
    """Count overlaps between two bed-like files."""
    intersect_file = Path("intersection.bed")
    intersect_args = [
        "-a",
        file_a,
        "-b",
        file_b,
        "-wa",
        "-u",
    ]
    return_code = (
        Cmd["bedtools"]["intersect"][intersect_args] > str(intersect_file)
    ) & RETCODE
    if return_code:
        error(f"Failed to intersect reads for {file_a} and {file_b}.")

    overlap_size = count_lines(intersect_file)

    intersect_file.unlink()
    return overlap_size


def rename_tagalign(file, name, tn5_shifted):
    """Rename tagAlign if it was transposon shifted."""
    if not tn5_shifted:
        return file
    tn5_tagalign = f"{name}_tn5.tagAlign.gz"
    Path(file).rename(tn5_tagalign)
    return tn5_tagalign


def process_narrow_peaks(file, cap_number):
    """Process narrow peak file for further analysis."""
    peaks = read_bed(
        bed_file=file,
        col_names=[
            "chromosome",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "signal_value",
            "p_value",
            "q_value",
            "peak",
        ],
    )

    peaks = peaks.sort_values(by=["p_value"], ascending=False)
    peaks["name"] = [f"Peak_{i}" for i in range(1, len(peaks.index) + 1)]
    peaks["start"] = peaks["start"].clip(lower=0)
    peaks["end"] = peaks["end"].clip(lower=0)

    if cap_number:
        peaks = peaks.head(cap_number)

    peaks.to_csv(file, sep="\t", index=False, header=False)


def shift_reads(tagalign, out_name, chromosome_sizes, frag_len, error):
    """Shift reads in a tagAlign and report the number of reads.

    Reads are shifted on both strands for half of the fragment length in
    opposite directions.
    """
    slop_args = [
        "-i",
        tagalign,
        "-g",
        chromosome_sizes,
        "-s",
        "-l",
        -int(frag_len / 2),
        "-r",
        int(frag_len / 2),
    ]
    return_code = (Cmd["bedtools"]["slop"][slop_args] > out_name) & RETCODE
    if return_code:
        error("Failed to shift reads in BED file.")

    shifted_df = read_bed(bed_file=out_name, col_names=BED_COLUMNS)
    shifted_df = shifted_df.loc[
        (shifted_df["start"] >= 0)
        & (shifted_df["end"] >= 0)
        & (shifted_df["start"] < shifted_df["end"])
    ]

    return len(shifted_df.index)


def get_frag_len(estimates):
    """Get the first fragment length estimate that is greater than 0."""
    for estimate in estimates:
        if estimate > 0:
            frag_len = estimate
            break
        else:
            frag_len = None
    return frag_len


class Macs2(Process):
    """Call peaks with MACS 2.0.

    Model-based Analysis of ChIP-Seq (MACS 2.0) is used to identify transcript
    factor binding sites using ChIP-Seq and related methods. MACS 2.0 captures
    the influence of genome complexity to evaluate the significance of enriched
    regions, and MACS improves the spatial resolution of binding sites through
    combining the information of both sequencing tag position and orientation.
    It has also an option to link nearby peaks together in order to call broad
    peaks. See [here](https://github.com/taoliu/MACS/) for more information.

    In addition to peak-calling, this process computes ChIP-Seq/CUT & RUN and
    ATAC-Seq QC metrics. Process returns a QC metrics report, fragment
    length estimation, and a deduplicated tagAlign file. QC report
    contains ENCODE 3 proposed QC metrics --
    [NRF](https://www.encodeproject.org/data-standards/terms/),
    [PBC bottlenecking coefficients, NSC, and RSC](https://genome.ucsc.edu/ENCODE/qualityMetrics.html#chipSeq).
    """

    slug = "macs2-callpeak"
    name = "MACS 2.0"
    process_type = "data:chipseq:callpeak:macs2"
    version = "4.8.2"
    category = "ChIP-seq"
    data_name = "{{ case|name|default('?') }}"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.CACHED
    entity = {"type": "sample", "input": "case"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/chipseq:6.0.0"}
        },
        "resources": {
            "cores": 4,
            "memory": 32768,
        },
    }

    class Input:
        """Input fields to process Macs2."""

        case = DataField(
            data_type="alignment:bam",
            label="Case (treatment)",
        )

        control = DataField(
            data_type="alignment:bam",
            label="Control (background)",
            required=False,
        )

        promoter = DataField(
            data_type="bed",
            label="Promoter regions BED file",
            required=False,
            description="BED file containing promoter regions (TSS+-1000bp "
            "for example). Needed to get the number of peaks and reads mapped "
            "to promoter regions.",
        )

        tagalign = BooleanField(
            label="Use tagAlign files",
            default=False,
            description="Use filtered tagAlign files as case "
            "(treatment) and control (background) samples. If extsize "
            "parameter is not set, run MACS using input's estimated fragment "
            "length.",
        )

        class PrepeakQC:
            """Pre-peak QC settings."""

            q_threshold = IntegerField(label="Quality filtering threshold", default=30)
            n_sub = IntegerField(label="Number of reads to subsample", default=15000000)
            tn5 = BooleanField(
                label="Tn5 shifting",
                default=False,
                description="Tn5 transposon shifting. Shift reads on '+' strand "
                "by 4bp and reads on '-' strand by 5bp.",
            )

            shift = IntegerField(
                label="User-defined cross-correlation peak strandshift",
                required=False,
                description="If defined, SPP tool will not try to estimate "
                "fragment length but will use the given value as "
                "fragment length.",
            )

        class Settings:
            """MACS2 settings."""

            format = StringField(
                label="Format of tag file",
                choices=[
                    ("BAM", "BAM"),
                    ("BAMPE", "BAMPE"),
                ],
                default="BAM",
                hidden="tagalign",
                description="This specifies the format of input files. For "
                "paired-end data the format dictates how MACS2 will treat "
                "mates. If the selected format is BAM, MACS2 will only keep "
                "the left mate (5' end) tag. However, when format BAMPE is "
                "selected, MACS2 will use actual insert sizes of pairs of "
                "reads to build fragment pileup, instead of building bimodal "
                "distribution plus and minus strand reads to predict fragment "
                "size.",
            )

            duplicates = StringField(
                label="Number of duplicates",
                choices=[
                    ("1", "1"),
                    ("auto", "auto"),
                    ("all", "all"),
                ],
                required=False,
                hidden="tagalign",
                description="It controls the MACS behavior towards duplicate "
                "tags at the exact same location -- the same coordination and "
                "the same strand. The 'auto' option makes MACS calculate the "
                "maximum tags at the exact same location based on binomial "
                "distribution using 1e-5 as pvalue cutoff and the 'all' "
                "option keeps all the tags. If an integer is given, at most "
                "this number of tags will be kept at the same location. The "
                "default is to keep one tag at the same location.",
            )

            duplicates_prepeak = StringField(
                label="Number of duplicates",
                choices=[
                    ("1", "1"),
                    ("auto", "auto"),
                    ("all", "all"),
                ],
                hidden="!tagalign",
                default="all",
                description="It controls the MACS behavior towards duplicate "
                "tags at the exact same location -- the same coordination and "
                "the same strand. The 'auto' option makes MACS calculate the "
                "maximum tags at the exact same location based on binomial "
                "distribution using 1e-5 as pvalue cutoff and the 'all' "
                "option keeps all the tags. If an integer is given, at most "
                "this number of tags will be kept at the same location. The "
                "default is to keep one tag at the same location.",
            )

            qvalue = FloatField(
                label="Q-value cutoff",
                required=False,
                disabled="settings.pvalue && settings.pvalue_prepeak",
                description="The q-value (minimum FDR) cutoff to call "
                "significant regions. Q-values are calculated from p-values "
                "using Benjamini-Hochberg procedure.",
            )

            pvalue = FloatField(
                label="P-value cutoff",
                disabled="settings.qvalue",
                hidden="tagalign",
                required=False,
                description="The p-value cutoff. If specified, MACS2 will use "
                "p-value instead of q-value cutoff.",
            )
            pvalue_prepeak = FloatField(
                label="P-value cutoff",
                default=0.00001,
                disabled="settings.qvalue",
                hidden="!tagalign || settings.qvalue",
                description="The p-value cutoff. If specified, MACS2 will use "
                "p-value instead of q-value cutoff.",
            )

            cap_num = IntegerField(
                label="Cap number of peaks by taking top N peaks",
                default=500000,
                disabled="settings.broad",
                description="To keep all peaks set value to 0.",
            )

            mfold_lower = IntegerField(
                label="MFOLD range (lower limit)",
                required=False,
                description="This parameter is used to select the regions "
                "within MFOLD range of high-confidence enrichment ratio "
                "against background to build model. The regions must be lower "
                "than upper limit, and higher than the lower limit of fold "
                "enrichment. DEFAULT:10,30 means using all regions not too "
                "low (>10) and not too high (<30) to build paired-peaks "
                "model. If MACS can not find more than 100 regions to build "
                "model, it will use the --extsize parameter to continue the "
                "peak detection ONLY if --fix-bimodal is set.",
            )

            mfold_upper = IntegerField(
                label="MFOLD range (upper limit)",
                required=False,
                description="This parameter is used to select the regions "
                "within MFOLD range of high-confidence enrichment ratio "
                "against background to build model. The regions must be lower "
                "than upper limit, and higher than the lower limit of fold "
                "enrichment. DEFAULT:10,30 means using all regions not too "
                "low (>10) and not too high (<30) to build paired-peaks "
                "model. If MACS can not find more than 100 regions to build "
                "model, it will use the --extsize parameter to continue the "
                "peak detection ONLY if --fix-bimodal is set.",
            )

            slocal = IntegerField(
                label="Small local region",
                required=False,
                description="Slocal and llocal parameters control which two "
                "levels of regions will be checked around the peak regions to "
                "calculate the maximum lambda as local lambda. By default, "
                "MACS considers 1000bp for small local region (--slocal), and "
                "10000bps for large local region (--llocal) which captures "
                "the bias from a long range effect like an open chromatin "
                "domain. You can tweak these according to your project. "
                "Remember that if the region is set too small, a sharp spike "
                "in the input data may kill the significant peak.",
            )

            llocal = IntegerField(
                label="Large local region",
                required=False,
                description="Slocal and llocal parameters control which two "
                "levels of regions will be checked around the peak regions to "
                "calculate the maximum lambda as local lambda. By default, "
                "MACS considers 1000bp for small local region (--slocal), and "
                "10000bps for large local region (--llocal) which captures "
                "the bias from a long range effect like an open chromatin "
                "domain. You can tweak these according to your project. "
                "Remember that if the region is set too small, a sharp spike "
                "in the input data may kill the significant peak.",
            )

            extsize = IntegerField(
                label="Extension size [--extsize]",
                required=False,
                description="While '--nomodel' is set, MACS uses this "
                "parameter to extend reads in 5'->3' direction to fix-sized "
                "fragments. For example, if the size of binding region for "
                "your transcription factor is 200 bp, and you want to bypass "
                "the model building by MACS, this parameter can be set as "
                "200. This option is only valid when --nomodel is set or when "
                "MACS fails to build model and --fix-bimodal is on.",
            )

            shift = IntegerField(
                label="Shift",
                required=False,
                hidden="settings.format == 'BAMPE'",
                description="Note, this is NOT the legacy --shiftsize option "
                "which is replaced by --extsize! You can set an arbitrary "
                "shift in bp here. Please Use discretion while setting it "
                "other than default value (0). When --nomodel is set, MACS "
                "will use this value to move cutting ends (5') then apply "
                "--extsize from 5' to 3' direction to extend them to "
                "fragments. When this value is negative, ends will be moved "
                "toward 3'->5' direction, otherwise 5'->3' direction. "
                "Recommended to keep it as default 0 for ChIP-Seq datasets, "
                "or -1 * half of EXTSIZE together with --extsize option for "
                "detecting enriched cutting loci such as certain DNAseI-Seq "
                "datasets. Note, you can't set values other than 0 if format "
                "is BAMPE for paired-end data. Default is 0.",
            )

            band_width = IntegerField(
                label="Band width",
                required=False,
                description="The band width which is used to scan the genome "
                "ONLY for model building. You can set this parameter as the "
                "sonication fragment size expected from wet experiment. The "
                "previous side effect on the peak detection process has been "
                "removed. So this parameter only affects the model building.",
            )

            nolambda = BooleanField(
                label="Use background lambda as local lambda",
                default=False,
                description="With this flag on, MACS will use the background "
                "lambda as local lambda. This means MACS will not consider "
                "the local bias at peak candidate regions.",
            )

            fix_bimodal = BooleanField(
                label="Turn on the auto paired-peak model process",
                default=False,
                description="Turn on the auto paired-peak model process. If "
                "it's set, when MACS failed to build paired model, it will "
                "use the nomodel settings, the '--extsize' parameter to "
                "extend each tag. If set, MACS will be terminated if "
                "paired-peak model has failed.",
            )

            nomodel = BooleanField(
                label="Bypass building the shifting model [--nomodel]",
                default=False,
                hidden="tagalign",
                description="While on, MACS will bypass building the shifting "
                "model.",
            )

            nomodel_prepeak = BooleanField(
                label="Bypass building the shifting model [--nomodel]",
                default=True,
                hidden="!tagalign",
                description="While on, MACS will bypass building the shifting "
                "model.",
            )

            down_sample = BooleanField(
                label="Down-sample",
                default=False,
                description="When set to true, random sampling method will "
                "scale down the bigger sample. By default, MACS uses linear "
                "scaling. This option will make the results unstable and "
                "unreproducible since each time, random reads would be "
                "selected, especially the numbers (pileup, pvalue, qvalue) "
                "would change.",
            )

            bedgraph = BooleanField(
                label="Save fragment pileup and control lambda",
                default=True,
                description="If this flag is on, MACS will store the fragment "
                "pileup, control lambda, -log10pvalue and -log10qvalue scores "
                "in bedGraph files. The bedGraph files will be stored in "
                "current directory named NAME+'_treat_pileup.bdg' for "
                "treatment data, NAME+'_control_lambda.bdg' for local lambda "
                "values from control, NAME+'_treat_pvalue.bdg' for Poisson "
                "pvalue scores (in -log10(pvalue) form), and "
                "NAME+'_treat_qvalue.bdg' for q-value scores from "
                "Benjamini-Hochberg-Yekutieli procedure.",
            )

            spmr = BooleanField(
                label="Save fragment pileup and control lambda",
                default=True,
                disabled="settings.bedgraph === false",
            )

            call_summits = BooleanField(
                label="Call summits [--call-summits]",
                default=False,
                description="MACS will now reanalyze the shape of signal "
                "profile (p or q-score depending on cutoff setting) to "
                "deconvolve subpeaks within each peak called from general "
                "procedure. It's highly recommended to detect adjacent "
                "binding events. While used, the output subpeaks of a big "
                "peak region will have the same peak boundaries, and "
                "different scores and peak summit positions.",
            )

            broad = BooleanField(
                label="Composite broad regions [--broad]",
                default=False,
                disabled="settings.call_summits === true",
                description="When this flag is on, MACS will try to composite "
                "broad regions in BED12 (a gene-model-like format) by "
                "putting nearby highly enriched regions into a broad region "
                "with loose cutoff. The broad region is controlled by another "
                "cutoff through --broad-cutoff. The maximum length of broad "
                "region length is 4 times of d from MACS.",
            )

            broad_cutoff = FloatField(
                label="Broad cutoff",
                required=False,
                disabled="settings.call_summits === true || settings.broad !== true",
                description="Cutoff for broad region. This option is not "
                "available unless --broad is set. If -p is set, this is a "
                "p-value cutoff, otherwise, it's a q-value cutoff. DEFAULT = "
                "0.1",
            )

        prepeakqc_settings = GroupField(PrepeakQC, label="Pre-peak QC settings")
        settings = GroupField(Settings, label="MACS2 settings")

    class Output:
        """Output fields to process Macs2."""

        called_peaks = FileField(label="Called peaks")
        narrow_peaks = FileField(label="Narrow peaks", required=False)
        chip_qc = FileField(label="QC report", required=False)
        case_prepeak_qc = FileField(label="Pre-peak QC report (case)")
        case_tagalign = FileField(label="Filtered tagAlign (case)")
        case_bam = FileField(label="Filtered BAM (case)")
        case_bai = FileField(label="Filtered BAM index (case)")
        control_prepeak_qc = FileField(
            label="Pre-peak QC report (control)", required=False
        )
        control_tagalign = FileField(
            label="Filtered tagAlign (control)", required=False
        )
        control_bam = FileField(label="Filtered BAM (control)", required=False)
        control_bai = FileField(label="Filtered BAM index (control)", required=False)
        narrow_peaks_bigbed_igv_ucsc = FileField(
            label="Narrow peaks (BigBed)", required=False
        )
        summits = FileField(label="Peak summits", required=False)
        summits_tbi_jbrowse = FileField(
            label="Peak summits tbi index for JBrowse", required=False
        )
        summits_bigbed_igv_ucsc = FileField(label="Summits (bigBed)", required=False)
        broad_peaks = FileField(label="Broad peaks", required=False)
        gappedPeak = FileField(label="Broad peaks (bed12/gappedPeak)", required=False)
        treat_pileup = FileField(label="Treatment pileup (bedGraph)", required=False)
        treat_pileup_bigwig = FileField(
            label="Treatment pileup (bigWig)", required=False
        )
        control_lambda = FileField(label="Control lambda (bedGraph)", required=False)
        control_lambda_bigwig = FileField(
            label="Control lambda (bigwig)", required=False
        )
        model = FileField(label="Model", required=False)
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""

        TMPDIR = os.environ.get("TMPDIR")

        # Allow Java to allocate 80% of the maximum memory available in
        # the container because the process itself uses some memory.
        java_memory = self.requirements.resources.memory // 1024 // 1.25
        os.environ["_JAVA_OPTIONS"] = f"-Xms256M -Xmx{java_memory}g"

        if inputs.settings.broad and inputs.settings.call_summits:
            self.error(
                "Composite broad regions [--broad] can't be combined with Call summits "
                "[--call-summits]."
            )

        if inputs.settings.shift and inputs.settings.format == "BAMPE":
            self.error(
                "Shift values other than 0 are not supported when the format is BAMPE."
            )

        if (
            inputs.settings.mfold_lower is not None
            and inputs.settings.mfold_upper is None
        ):
            self.error(
                "MFOLD range should be set both for upper and lower limit, but only the lower "
                "limit was set."
            )

        if (
            inputs.settings.mfold_lower is None
            and inputs.settings.mfold_upper is not None
        ):
            self.error(
                "MFOLD range should be set both for upper and lower limit, but only the upper "
                "limit was set."
            )

        if inputs.control:
            if inputs.case.output.species != inputs.control.output.species:
                self.error(
                    "All input BAM files must share the same genome species information. BAM "
                    f"{inputs.case.name} has {inputs.case.output.species} while "
                    f"{inputs.control.name} has {inputs.control.output.species} species "
                    "information."
                )
            if inputs.case.output.build != inputs.control.output.build:
                self.error(
                    "All input BAM files must share the same genome build information. BAM"
                    f"{inputs.case.name} has {inputs.case.output.build} while "
                    f"{inputs.control.name} has {inputs.control.output.build} build information."
                )

        if inputs.promoter:
            if inputs.promoter.output.species != inputs.case.output.species:
                self.error(
                    "All input files must share the same genome species information. Case BAM "
                    f"has {inputs.case.output.species} while promoter BED has "
                    f"{inputs.promoter.output.species} species information."
                )
            if inputs.promoter.output.build != inputs.case.output.build:
                self.error(
                    "All input files must share the same genome build information. Case BAM "
                    f"has {inputs.case.output.build} while promoter BED has"
                    f"{inputs.promoter.output.build} build information."
                )

        try:
            gsize = SPECIES_GSIZES[inputs.case.output.species]
        except KeyError:
            self.error(
                f"Species {inputs.case.output.species} is not supported by the MACS 2.0 process. "
                f"Supported species are: {', '.join(SPECIES_GSIZES.keys())}"
            )

        self.progress(0.1)

        samples = [inputs.case, inputs.control] if inputs.control else [inputs.case]
        for alignment in samples:
            bam_path = Path(alignment.output.bam.path)
            is_control = alignment == inputs.control if inputs.control else False
            is_paired = test_paired_end(bam_file=bam_path, error=self.error)
            if not is_paired and inputs.settings.format == "BAMPE":
                self.error(
                    "No paired-end reads were detected but BAMPE format was selected."
                )

            name = f"{bam_path.stem}_background" if is_control else f"{bam_path.stem}"
            # Set case name and layout for post-qc analysis.
            if not is_control:
                case_name = name
                case_paired = is_paired

            if alignment.type.startswith("data:alignment:bam:bwaaln"):
                with open(alignment.output.stats.path, "r") as handle:
                    report = handle.read()
                metrics = parse_flagstat(report=report, warning=self.warning)
            elif alignment.type.startswith("data:alignment:bam:bowtie2"):
                return_code, n_mapped, _ = Cmd["samtools"][
                    "view", "-c", str(bam_path)
                ] & TEE(retcode=None)
                if return_code:
                    n_mapped = "NA"
                    self.warning(
                        "Failed to determine number of mapped reads based on Bowtie2 alignment."
                    )
                else:
                    n_mapped = int(n_mapped.strip())

                metrics = parse_bowtie2_report(
                    report_path=alignment.output.stats.path,
                    n_mapped=n_mapped,
                    warning=self.warning,
                )
            else:
                return_code, report, _ = Cmd["samtools"]["flagstat"][
                    str(bam_path)
                ] & TEE(retcode=None)
                if return_code:
                    self.error(f"Samtools flagstat failed for {name}")
                metrics = parse_flagstat(report=report, warning=self.warning)

            # Remove unmapped reads, not primary alignments, and reads below minimum mapping quality
            filtered_bam = f"{name}_filtered.bam"
            filter_bam(
                bam_file=str(bam_path),
                out_bam=filtered_bam,
                min_quality=inputs.prepeakqc_settings.q_threshold,
                is_paired=is_paired,
                name=name,
                error=self.error,
            )

            if is_paired:
                temp_bam = f"{name}_tmp.bam"
                name_sort_bam(in_bam=filtered_bam, out_bam=temp_bam, error=self.error)

                fixmate_bam = f"{case_name}_fixmate.bam"
                fixmate_params = ["-r", temp_bam, fixmate_bam]
                return_code, stdout, stderr = Cmd["samtools"]["fixmate"][
                    fixmate_params
                ] & TEE(retcode=None)

                if return_code:
                    print(stdout, stderr)
                    self.error(f"Samtools fixmate failed for {alignment.name}.")

                filter_bam(
                    bam_file=fixmate_bam,
                    out_bam=temp_bam,
                    min_quality=0,
                    is_paired=is_paired,
                    name=name,
                    error=self.error,
                )

                # Sort by position
                sort_params = ["-o", filtered_bam, temp_bam]
                return_code, stdout, stderr = Cmd["samtools"]["sort"][
                    sort_params
                ] & TEE(retcode=None)

                if return_code:
                    print(stdout, stderr)
                    self.error(f"Samtools sorting failed for {alignment.name}.")

            markdup_bam = f"{name}_markdup.bam"
            markdup_stats = f"{name}_duplicates_metrics.txt"
            markdup_params = [
                f"INPUT={filtered_bam}",
                f"OUTPUT={markdup_bam}",
                f"METRICS_FILE={markdup_stats}",
                "VALIDATION_STRINGENCY=LENIENT",
                "ASSUME_SORTED=true",
                "REMOVE_DUPLICATES=false",
                f"TMP_DIR={TMPDIR}",
            ]

            return_code, stdout, stderr = Cmd["java"]["-jar"][
                "/opt/broadinstitute/picard-tools/picard.jar"
            ]["MarkDuplicates"][markdup_params] & TEE(retcode=None)

            if return_code:
                print(stdout, stderr)
                self.error(
                    f"Picard-tools MarkDuplicates processing failed for {alignment.name}."
                )

            metrics = merge_dict(
                metrics,
                parse_markdup(markdup_stats),
            )

            # Name sort for PBC bottlenecking calculation
            if is_paired:
                name_sort_bam(in_bam=markdup_bam, out_bam=temp_bam, error=self.error)

            alignment_bed = f"{name}.bed"
            metrics = merge_dict(
                metrics,
                get_pbc_metrics(
                    bam_file=temp_bam if is_paired else markdup_bam,
                    name=name,
                    out_bed=alignment_bed,
                    is_paired=is_paired,
                    error=self.error,
                ),
            )

            # Remove unmapped reads, not primary alignments, and duplicate reads
            filter_bam(
                bam_file=markdup_bam,
                out_bam=filtered_bam,
                min_quality=0,
                is_paired=is_paired,
                name=name,
                error=self.error,
            )

            # Remove redundant bam with marked duplicates.
            Path(markdup_bam).unlink()

            if is_paired:
                name_sort_bam(in_bam=filtered_bam, out_bam=temp_bam, error=self.error)

            convert_bam(
                bam_file=temp_bam if is_paired else filtered_bam,
                name=name,
                out_bed=alignment_bed,
                is_paired=is_paired,
                error=self.error,
                first_mate=True,
            )

            # Remove redundant temporary bam.
            if is_paired:
                Path(temp_bam).unlink()

            tagalign = f"{name}.tagAlign.gz"
            prepare_tagalign(
                bed_file=alignment_bed,
                tagalign_file=tagalign,
                is_paired=is_paired,
                tn5_shift=inputs.prepeakqc_settings.tn5,
                compression="gzip",
            )

            subsampled_tagalign = f"{name}_subsampled.tagAlign.gz"
            subsample_tagalign(
                tagalign_file=tagalign,
                out_file=subsampled_tagalign,
                is_paired=is_paired,
                n_sub=inputs.prepeakqc_settings.n_sub,
                compression="gzip",
            )

            cross_correlation_report = f"{name}_cc_score.txt"
            spp_params = [
                f"-c={subsampled_tagalign}",
                f"-p={self.requirements.resources.cores}",
                "-filtchr=chrM",
                f"-out={cross_correlation_report}",
            ]

            if inputs.prepeakqc_settings.shift is not None:
                spp_params.append(f"-speak={inputs.prepeakqc_settings.shift}")

            return_code, stdout, stderr = Cmd["spp"][spp_params] & TEE(retcode=None)

            if return_code:
                print(stdout, stderr)
                self.error(f"SPP processing failed for {alignment.name}.")

            cc_report = pd.read_csv(
                cross_correlation_report,
                sep="\t",
                index_col=0,
                names=SPP_HEADER,
                dtype="str",
            )

            # Some columns have top 3 predictions and only the first one is needed.
            cc_metrics = cc_report.apply(lambda col: col.str.split(",").str[0])

            # Use case's estimated fragment length only.
            if not is_control:
                fraglen_estimates = cc_report.loc[
                    subsampled_tagalign, "Est. Fragment Len."
                ]

                estimate_list = [int(x) for x in fraglen_estimates.split(",")]
                frag_len = get_frag_len(estimates=estimate_list)

                if inputs.tagalign and inputs.settings.extsize is None:
                    if not frag_len:
                        self.error(
                            "Failed to estimate fragment length. No estimates were larger than "
                            f"zero. The top estimates were: {fraglen_estimates}. Please manually "
                            "specify the Extension size [--extsize] parameter."
                        )
                    elif frag_len != estimate_list[0]:
                        self.warning(
                            "SPP estimated negative fragment length which can not be used by "
                            f"MACS2. Using {frag_len} from the top estimates "
                            f"({fraglen_estimates}) as the estimate of extension size [--extsize] "
                            "for MACS2."
                        )

                # When not using tagAlign as an input we need to keep the first
                # estimate for post-peak QC steps.
                if (
                    not inputs.tagalign
                    and inputs.settings.extsize is None
                    and frag_len is None
                ):
                    frag_len = estimate_list[0]
                    self.warning(
                        "No fragment length estimate was greater than 0. Using the first "
                        f"estimate: {frag_len} for read shifting in post-peak QC."
                    )

            metrics = merge_dict(metrics, cc_metrics.loc[subsampled_tagalign].to_dict())

            prepeak_qc_report = f"{name}_prepeak_qc_report.txt"
            save_prepeak_qc(metrics=metrics, out_file=prepeak_qc_report)

            filtered_bai = f"{filtered_bam}.bai"
            return_code, stdout, stderr = Cmd["samtools"]["index"][
                filtered_bam, filtered_bai
            ] & TEE(retcode=None)

            if return_code:
                print(stdout, stderr)
                self.error(f"Samtools index failed for {filtered_bam}.")

            if is_control:
                outputs.control_prepeak_qc = prepeak_qc_report
                outputs.control_bam = filtered_bam
                outputs.control_bai = filtered_bai

                control_tagalign = rename_tagalign(
                    file=tagalign, name=name, tn5_shifted=inputs.prepeakqc_settings.tn5
                )
                outputs.control_tagalign = control_tagalign

            else:
                outputs.case_prepeak_qc = prepeak_qc_report
                outputs.case_bam = filtered_bam
                outputs.case_bai = filtered_bai

                case_tagalign = rename_tagalign(
                    file=tagalign, name=name, tn5_shifted=inputs.prepeakqc_settings.tn5
                )
                outputs.case_tagalign = case_tagalign

        self.progress(0.6)

        callpeak_params = [
            "-t",
            case_tagalign if inputs.tagalign else inputs.case.output.bam.path,
            "-f",
            "BED" if inputs.tagalign else inputs.settings.format,
            "-n",
            case_name,
            "--gsize",
            gsize,
            "--verbose",
            3,
            "--outdir",
            "./",
        ]

        if inputs.control:
            callpeak_params.extend(
                [
                    "-c",
                    (
                        control_tagalign
                        if inputs.tagalign
                        else inputs.control.output.bam.path
                    ),
                ]
            )

        if inputs.settings.duplicates or inputs.tagalign:
            callpeak_params.extend(
                [
                    "--keep-dup",
                    inputs.settings.duplicates or inputs.settings.duplicates_prepeak,
                ]
            )

        if inputs.settings.qvalue is not None and inputs.settings.pvalue is None:
            callpeak_params.extend(["-q", inputs.settings.qvalue])

        if inputs.settings.pvalue is not None and inputs.settings.qvalue is None:
            callpeak_params.extend(["-p", inputs.settings.pvalue])
        elif inputs.tagalign and inputs.settings.qvalue is None:
            callpeak_params.extend(["-p", inputs.settings.pvalue_prepeak])

        if inputs.settings.mfold_lower is not None:
            callpeak_params.extend(
                [
                    "-m",
                    inputs.settings.mfold_lower,
                    inputs.settings.mfold_upper,
                ]
            )

        if inputs.settings.nolambda:
            callpeak_params.append("--nolambda")

        if inputs.settings.slocal is not None:
            callpeak_params.extend(["--slocal", inputs.settings.slocal])

        if inputs.settings.llocal is not None:
            callpeak_params.extend(["--llocal", inputs.settings.llocal])

        if inputs.settings.fix_bimodal:
            callpeak_params.append("--fix-bimodal")

        if inputs.settings.nomodel or (
            inputs.settings.nomodel_prepeak and inputs.tagalign
        ):
            callpeak_params.append("--nomodel")

        if inputs.settings.extsize is not None:
            callpeak_params.extend(["--extsize", inputs.settings.extsize])
        elif inputs.tagalign:
            callpeak_params.extend(["--extsize", frag_len])

        if inputs.settings.shift is not None:
            callpeak_params.extend(["--shift", inputs.settings.shift])

        if inputs.settings.band_width is not None:
            callpeak_params.extend(["--bw", inputs.settings.band_width])

        if inputs.settings.broad:
            callpeak_params.append("--broad")
            if inputs.settings.broad_cutoff is not None:
                callpeak_params.extend(["--broad-cutoff", inputs.settings.broad_cutoff])

        if inputs.settings.down_sample:
            callpeak_params.append("--down-sample")

        if inputs.settings.bedgraph:
            callpeak_params.append("-B")
            if inputs.settings.spmr:
                callpeak_params.append("--SPMR")

        if inputs.settings.call_summits:
            callpeak_params.append("--call-summits")

        return_code, stdout, stderr = Cmd["macs2"]["callpeak"][callpeak_params] & TEE(
            retcode=None
        )

        if return_code:
            print(stdout, stderr)
            self.error("MACS2 processing failed.")

        print(stderr)
        warning_lines = [line for line in stderr.split("\n") if "WARNING" in line]
        if warning_lines:
            macs2_warning = ", ".join(warning_lines[:3]) + (
                ", ..." if len(warning_lines) > 3 else ""
            )
            self.warning(f"MACS2 reported warnings: {macs2_warning}")

        self.progress(0.7)

        model_script = f"{case_name}_model.r"
        if Path(model_script).is_file():
            return_code, _, _ = Cmd["Rscript"][model_script] & TEE(retcode=None)

            if return_code:
                self.error(f"Running R script {model_script} failed.")

            outputs.model = f"{case_name}_model.pdf"

        self.progress(0.8)

        outputs.called_peaks = f"{case_name}_peaks.xls"

        # Get chromosome sizes file for bed to BigBed transformation.
        return_code, idxstats, stderr = Cmd["samtools"]["idxstats"][
            inputs.case.output.bam.path
        ] & TEE(retcode=None)

        if return_code:
            print(stderr)
            self.error("Samtools idxstats failed")

        chromosome_sizes = "chrom.sizes"
        with open(chromosome_sizes, "w") as sizes_file:
            sizes_file.write("\n".join(idxstats.split("\n")[:-1]))

        if inputs.settings.broad:
            outputs.broad_peaks = f"{case_name}_peaks.broadPeak"
            outputs.gappedPeak = f"{case_name}_peaks.gappedPeak"
        else:
            # Maximize 5th column of narrowPeak and summits files to 1000.
            narrow_peak = f"{case_name}_peaks.narrowPeak"
            corrected_narrow_peak = f"corrected_{narrow_peak}"
            correct_bed_file(
                in_bed=narrow_peak, out_bed=corrected_narrow_peak, error=self.error
            )

            summits = f"{case_name}_summits.bed"
            corrected_summits = f"corrected_{summits}"
            correct_bed_file(
                in_bed=summits, out_bed=corrected_summits, error=self.error
            )

            narrow_peak_bb = f"{case_name}_peaks_narrowPeak.bb"
            create_big_bed(
                in_file=corrected_narrow_peak,
                out_bb=narrow_peak_bb,
                chromosome_sizes=chromosome_sizes,
                file_type="narrowPeak",
                error=self.error,
            )

            summits_bb = f"{case_name}_summits.bb"
            create_big_bed(
                in_file=corrected_summits,
                out_bb=summits_bb,
                chromosome_sizes=chromosome_sizes,
                file_type="summits.bed",
                error=self.error,
            )

            # Create tabix index for summits.bed file for JBrowse.
            summits_gz = f"{summits}.gz"
            (Cmd["bgzip"]["-c", summits] > summits_gz)()

            return_code, stdout, stderr = Cmd["tabix"]["-p", "bed", summits_gz] & TEE(
                retcode=None
            )

            if return_code:
                print(stdout, stderr)
                self.error("Summits.bed tabix processing for JBrowse failed.")

            outputs.narrow_peaks_bigbed_igv_ucsc = narrow_peak_bb
            outputs.summits = summits_gz
            outputs.summits_tbi_jbrowse = f"{summits_gz}.tbi"
            outputs.summits_bigbed_igv_ucsc = summits_bb

        if inputs.settings.bedgraph:
            pileup_bdg = f"{case_name}_treat_pileup.bdg"
            outputs.treat_pileup = pileup_bdg

            pileup_bw = f"{case_name}_treat_pileup.bw"
            create_bigwig(
                in_bdg=pileup_bdg,
                out_bw=pileup_bw,
                chromosome_sizes=chromosome_sizes,
                file_type="treat_pileup.bgd bedGraph",
                error=self.error,
            )

            outputs.treat_pileup_bigwig = pileup_bw

            control_lambda_bdg = f"{case_name}_control_lambda.bdg"
            outputs.control_lambda = control_lambda_bdg

            control_lambda_bw = f"{case_name}_control_lambda.bw"
            create_bigwig(
                in_bdg=control_lambda_bdg,
                out_bw=control_lambda_bw,
                chromosome_sizes=chromosome_sizes,
                file_type="control_lambda.bgd bedGraph",
                error=self.error,
            )

            outputs.control_lambda_bigwig = control_lambda_bw

        if not inputs.settings.broad:
            process_narrow_peaks(file=narrow_peak, cap_number=inputs.settings.cap_num)

            # Use either filtered tagAlign file or compute regular tagAlign file for the case sample.
            if not inputs.tagalign:
                if case_paired:
                    temp_case_bam = f"{case_name}_tmp.bam"
                    name_sort_bam(
                        in_bam=inputs.case.output.bam.path,
                        out_bam=temp_case_bam,
                        error=self.error,
                    )

                    fixmate_bam = f"{case_name}_fixmate.bam"
                    fixmate_params = ["-r", temp_case_bam, fixmate_bam]
                    return_code, stdout, stderr = Cmd["samtools"]["fixmate"][
                        fixmate_params
                    ] & TEE(retcode=None)

                    if return_code:
                        print(stdout, stderr)
                        self.error(f"Samtools fixmate failed for {case_name}.")

                    args = ["-f", "2", "-o", temp_case_bam, fixmate_bam]
                    return_code, stdout, stderr = Cmd["samtools"]["view"][args] & TEE(
                        retcode=None
                    )
                    if return_code:
                        print(stdout, stderr)
                        self.error(f"Samtools view failed for {case_name}.")

                temp_case_bed = f"{case_name}_tmp.bed"
                bam_to_convert = (
                    temp_case_bam if is_paired else inputs.case.output.bam.path
                )

                convert_bam(
                    bam_file=bam_to_convert,
                    out_bed=temp_case_bed,
                    name=case_name,
                    is_paired=case_paired,
                    error=self.error,
                )

                case_tagalign = f"{case_name}_nonfiltered.tagAlign"
                prepare_tagalign(
                    bed_file=temp_case_bed,
                    tagalign_file=case_tagalign,
                    is_paired=is_paired,
                    tn5_shift=False,
                )

            # Shift reads on both strands for half of the fragment length in opposite directions.
            shifted_tagaling = f"{case_name}_shifted.tagAlign"
            n_reads = shift_reads(
                tagalign=case_tagalign,
                out_name=shifted_tagaling,
                chromosome_sizes=chromosome_sizes,
                frag_len=inputs.settings.extsize or frag_len,
                error=self.error,
            )

            # Calculate post-peakcall QC metrics
            post_peak_qc = {
                "FRiP": "N/A",
                "NUMBER_OF_PEAKS": count_lines(narrow_peak),
                "NUMBER_OF_READS_IN_PROMOTERS": "N/A",
                "FRACTION_OF_READS_IN_PROMOTERS": "N/A",
                "NUMBER_OF_PEAKS_IN_PROMOTERS": "N/A",
                "FRACTION_OF_PEAKS_IN_PROMOTERS": "N/A",
            }

            n_reads_peaks = count_overlap(
                file_a=shifted_tagaling, file_b=narrow_peak, error=self.error
            )

            post_peak_qc["FRiP"] = round(n_reads_peaks / n_reads, 3)
            if inputs.promoter:
                post_peak_qc["NUMBER_OF_READS_IN_PROMOTERS"] = count_overlap(
                    file_a=shifted_tagaling,
                    file_b=inputs.promoter.output.bed.path,
                    error=self.error,
                )
                post_peak_qc["FRACTION_OF_READS_IN_PROMOTERS"] = round(
                    post_peak_qc["NUMBER_OF_READS_IN_PROMOTERS"] / n_reads, 3
                )

                post_peak_qc["NUMBER_OF_PEAKS_IN_PROMOTERS"] = count_overlap(
                    file_a=narrow_peak,
                    file_b=inputs.promoter.output.bed.path,
                    error=self.error,
                )
                post_peak_qc["FRACTION_OF_PEAKS_IN_PROMOTERS"] = round(
                    post_peak_qc["NUMBER_OF_PEAKS_IN_PROMOTERS"]
                    / post_peak_qc["NUMBER_OF_PEAKS"],
                    3,
                )

            qc_file = f"{case_name}_postpeak_qc_report.txt"
            pd.Series(post_peak_qc, dtype="object").to_frame().T.to_csv(
                qc_file, sep="\t", index=False
            )

            return_code, _, _ = Cmd["gzip"][narrow_peak] & TEE(retcode=None)
            if return_code:
                self.error("Compression of narrowPeaks file failed.")

            outputs.narrow_peaks = f"{narrow_peak}.gz"
            outputs.chip_qc = qc_file

        outputs.build = inputs.case.output.build
        outputs.species = inputs.case.output.species
