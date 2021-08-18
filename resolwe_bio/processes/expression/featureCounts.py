"""featureCounts."""
import gzip
import io
import json
import shutil
from pathlib import Path

import pandas as pd
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    FloatField,
    GroupField,
    IntegerField,
    JsonField,
    SchedulingClass,
    StringField,
)

from resolwe_bio.process.runtime import ProcessBio


def get_gene_counts(infile, outfile, sample_name):
    """Rename columns and extract gene counts from featureCounts output file."""
    exp = pd.read_csv(
        infile,
        sep="\t",
        skip_blank_lines=True,
        header=1,
        usecols=["Geneid", sample_name],
        index_col="Geneid",
        dtype={
            "Geneid": str,
            sample_name: int,
        },
        squeeze=True,
    )
    return exp.to_csv(
        outfile,
        index_label="FEATURE_ID",
        header=["EXPRESSION"],
        sep="\t",
    )


def get_gene_lenghts(infile, outfile):
    """Rename columns and extract feature lengths from featureCounts output file."""
    exp = pd.read_csv(
        infile,
        sep="\t",
        skip_blank_lines=True,
        header=1,
        usecols=["Geneid", "Length"],
        index_col="Geneid",
        dtype={
            "Geneid": str,
            "Length": int,
        },
        squeeze=True,
    )
    return exp.to_csv(
        outfile,
        index_label="FEATURE_ID",
        header=["GENE_LENGTHS"],
        sep="\t",
    )


def rename_columns_and_compress(infile, outfile):
    """Rename columns and compress the expression files."""
    exp = pd.read_csv(
        infile,
        sep="\t",
        usecols=["FEATURE_ID", "EXPRESSION"],
        index_col="FEATURE_ID",
        squeeze=True,
        float_precision="round_trip",
    )
    return exp.to_csv(
        outfile, index_label="Gene", header=["Expression"], sep="\t", compression="gzip"
    )


def exp_to_df(exp_file, exp_type):
    """Prepare expression file for gene sets merging."""
    with open(exp_file) as exp:
        df = pd.read_csv(exp, sep="\t", float_precision="round_trip")
        df.rename(
            index=str,
            columns={
                "EXPRESSION": exp_type,
            },
            inplace=True,
        )
        # Cast FEATURE_ID column to string
        df["FEATURE_ID"] = df["FEATURE_ID"].astype("str")

    return df


def prepare_expression_set(rc, tpm, cpm, feature_dict, outfile_name):
    """Prepare expression set output data."""
    rc_exp = exp_to_df(rc, "RAW_COUNT")
    tpm_exp = exp_to_df(tpm, "TPM")
    cpm_exp = exp_to_df(cpm, "CPM")
    rc_exp["GENE_SYMBOL"] = rc_exp["FEATURE_ID"].map(feature_dict)
    input_features = rc_exp["FEATURE_ID"].tolist()
    # Check if all of the input feature IDs could be mapped to the gene symbols
    if not all(f_id in feature_dict for f_id in input_features):
        print(
            f"{sum(rc_exp.isnull().values.ravel())} feature(s) "
            f"could not be mapped to the associated feature symbols."
        )
    # Merge with normalized expression values
    exp_set = rc_exp.merge(tpm_exp, on="FEATURE_ID")
    exp_set = exp_set.merge(cpm_exp, on="FEATURE_ID")
    # Reorder columns
    columns = ["FEATURE_ID", "GENE_SYMBOL", "RAW_COUNT", "TPM", "CPM"]
    exp_set = exp_set[columns]
    # Replace NaN values with empty string
    exp_set.fillna("", inplace=True)

    # Write to file
    exp_set.to_csv(
        outfile_name + ".txt.gz",
        header=True,
        index=False,
        sep="\t",
        compression="gzip",
    )

    # Write to JSON
    df_dict = exp_set.set_index("FEATURE_ID").to_dict(orient="index")
    with open(outfile_name + ".json", "w") as f:
        json.dump({"genes": df_dict}, f, allow_nan=False)


def expression_to_storage(rc_input, rc_output):
    """Convert expressions file to JSON format."""

    def isfloat(value):
        """Check if value is float."""
        try:
            float(value)
            return True
        except ValueError:
            return False

    with io.TextIOWrapper(io.BufferedReader(gzip.open(rc_input))) as f:
        # Split lines by tabs
        # Ignore lines without a number in second column
        # Build a dictionary of gene-expression pairs
        exp = {
            "genes": {
                gene_exp[0]: float(gene_exp[1])
                for gene_exp in (l.split("\t") for l in f)
                if len(gene_exp) == 2 and isfloat(gene_exp[1])
            }
        }

    with open(file=rc_output, mode="wt") as f:
        json.dump(exp, f)

    return rc_output


class featureCounts(ProcessBio):
    """Quantify sequencing reads aligned to genomic features.

    featureCounts is a highly efficient general-purpose read summarization
    program that counts aligned reads on genomic features such as genes, exons,
    promoter, gene bodies, genomic bins and chromosomal locations. It can be
    used to count both RNA-seq and genomic DNA-seq reads. See the
    [official website](http://bioinf.wehi.edu.au/featureCounts/) and the
    [introductory paper](https://academic.oup.com/bioinformatics/article/30/7/923/232889)
    for more information.

    featureCounts output includes raw counts and normalized (TPM/CPM) expression values.
    Normalized expression values are computed using rnanorm Python package under
    union-exon gene length model.
    """

    slug = "feature_counts"
    name = "featureCounts"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/s4q6j6e8/resolwebio/rnaseq:5.11.0",
            },
        },
        "resources": {
            "cores": 4,
            "memory": 32768,
            "network": True,
        },
    }
    data_name = "{{ aligned_reads|sample_name|default('?') }}"
    version = "5.0.0"
    process_type = "data:expression:featurecounts"
    category = "Quantify"
    entity = {
        "type": "sample",
    }
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        aligned_reads = DataField("alignment:bam", label="Aligned reads")

        annotation = DataField(
            "annotation",
            label="Annotation",
            description="GTF and GFF3 annotation formats are supported",
        )

        feature_class = StringField(
            label="Feature class",
            default="exon",
            description="Feature class (3rd column in GTF/GFF3 file) to be used. All other features will be ignored.",
        )

        id_attribute = StringField(
            label="ID attribute",
            allow_custom_choice=True,
            default="gene_id",
            choices=[
                ("gene_id", "gene_id"),
                ("transcript_id", "transcript_id"),
                ("ID", "ID"),
                ("geneid", "geneid"),
                ("Name", "Name"),
            ],
            description="GTF/GFF3 attribute to be used as feature ID. Several GTF/GFF3 lines "
            "with the same feature ID will be considered as parts of the same "
            "feature. The feature ID is used to identify the counts in the "
            "output table. In GTF files this is usually 'gene_id', in GFF3 files "
            "this is often 'ID', and 'transcript_id' is frequently a valid "
            "choice for both annotation formats.",
        )

        feature_type = StringField(
            label="Feature type",
            default="gene",
            choices=[
                ("gene", "gene"),
                ("transcript", "transcript"),
            ],
            description="The type of feature the quantification program summarizes over "
            "(e.g. gene or transcript-level analysis). The value of this "
            "parameter needs to be chosen in line with ID attribute input choice.",
        )

        normalization_type = StringField(
            label="Normalization type",
            default="TPM",
            choices=[
                ("TPM", "TPM"),
                ("CPM", "CPM"),
            ],
            description="The default expression normalization type.",
        )

        assay_type = StringField(
            label="Assay type",
            default="non_specific",
            choices=[
                ("non_specific", "Strand non-specific"),
                ("forward", "Strand-specific forward"),
                ("reverse", "Strand-specific reverse"),
                ("auto", "Detect automatically"),
            ],
            description="Indicate if strand-specific read counting should be performed. "
            "For paired-end reads, strand of the first read is taken as the strand "
            "of the whole fragment. FLAG field is used to tell if a read is "
            "first or second read in a pair. Automated strand detection is enabled "
            "using the [Salmon](https://salmon.readthedocs.io/en/latest/library_type.html) "
            "tool's build-in functionality. To use this option, cDNA (transcriptome) "
            "index file crated using the Salmon indexing tool must be provided",
        )

        cdna_index = DataField(
            "index:salmon",
            label="Salmon index file",
            required=False,
            hidden="assay_type != 'auto'",
            description="Transcriptome index file created using the Salmon indexing tool. "
            "cDNA (transcriptome) sequences used for index file creation must be "
            "derived from the same species as the input sequencing reads to "
            "obtain the reliable analysis results.",
        )

        n_reads = IntegerField(
            label="Number of reads in subsampled alignment file",
            default=5000000,
            hidden="assay_type != 'auto'",
            description="Alignment (.bam) file subsample size to detect "
            "strandedness. Increase the number of reads make automatic "
            "detection more reliable. Decrease the number of reads to "
            "make automatic detection run faster.",
        )

        advanced = BooleanField(
            label="Show advanced options",
            description="Inspect and modify parameters.",
            default=False,
        )

        class General:
            """General options."""

            count_features = BooleanField(
                label="Perform read counting at feature level",
                description="Count reads for exons rather than genes.",
                default=False,
            )

            by_read_group = BooleanField(
                label="Assign reads by read group",
                description="RG tag is required to be present in the input BAM files.",
                default=False,
            )

            count_long_reads = BooleanField(
                label="Count long reads such as Nanopore and PacBio reads",
                default=False,
            )

            count_multi_mapping_reads = BooleanField(
                label="Count multi-mapping reads",
                description="For a multi-mapping read, all its reported alignments will be "
                "counted. The 'NH' tag in BAM input is used to detect multi-mapping reads.",
                default=False,
            )

            fraction = BooleanField(
                label="Assign fractional counts to features",
                description="This option must be used together with 'Count multi-mapping "
                "reads' or 'Assign reads to all their overlapping features or "
                "meta-features' or both. When 'Count multi-mapping reads' is "
                "checked, each reported alignment from a multi-mapping read "
                "(identified via 'NH' tag) will carry a count of 1 / x, instead "
                "of 1 (one), where x is the total number of alignments reported "
                "for the same read. When 'Assign reads to all their overlapping "
                "features or meta-features' is checked, each overlapping "
                "feature will receive a count of 1 / y, where y is the total "
                "number of features overlapping with the read. When both 'Count "
                "multi-mapping reads' and 'Assign reads to all their overlapping "
                "features or meta-features' are specified, each alignment will "
                "carry a count of 1 / (x * y).",
                default=False,
                disabled="!(general.count_multi_mapping_reads || overlap.allow_multi_overlap)",
            )

        class Overlap:
            """Overlap options."""

            allow_multi_overlap = BooleanField(
                label="Assign reads to all their overlapping features or meta-features",
                default=False,
            )

            min_overlap = IntegerField(
                label="Minimum number of overlapping bases in a read that is required for read assignment",
                default=1,
                description="Number of overlapping bases is counted from both reads if "
                "paired-end. If a negative value is provided, then a gap of up "
                "to specified size will be allowed between read and the feature "
                "that the read is assigned to.",
            )

            frac_overlap = FloatField(
                label="Minimum fraction of overlapping bases in a read that is required for read assignment",
                default=0.0,
                description="Value should be within range [0, 1]. Number of overlapping bases "
                "is counted from both reads if paired end. Both this "
                "option and 'Minimum number of overlapping bases in a read "
                "that is required for read assignment' need to be satisfied "
                "for read assignment.",
            )

            frac_overlap_feature = FloatField(
                label="Minimum fraction of overlapping bases included in a feature that is "
                "required for overlapping with a read or a read pair",
                default=0.0,
                description="Value should be within range [0, 1].",
            )

            largest_overlap = BooleanField(
                label="Assign reads to a feature or meta-feature that has the largest "
                "number of overlapping bases",
                default=False,
            )

            read_extension_5 = IntegerField(
                label="Number of bases to extend reads upstream by from their 5' end",
                default=0,
            )

            read_extension_3 = IntegerField(
                label="Number of bases to extend reads upstream by from their 3' end",
                default=0,
            )

            read_to_pos = IntegerField(
                label="Reduce reads to their 5'-most or 3'-most base",
                required=False,
                description="Read counting is performed based on the single base the read "
                "is reduced to.",
            )

        class ReadFiltering:
            """Read filtering."""

            min_mqs = IntegerField(
                label="Minimum mapping quality score",
                default=0,
                description="The minimum mapping quality score a read must satisfy in order "
                "to be counted. For paired-end reads, at least one end should "
                "satisfy this criterion.",
            )

            split_only = BooleanField(
                label="Count only split alignments",
                default=False,
            )

            non_split_only = BooleanField(
                label="Count only non-split alignments",
                default=False,
            )

            primary = BooleanField(
                label="Count only primary alignments",
                default=False,
                description="Primary alignments are identified using bit 0x100 in BAM "
                "FLAG field.",
            )

            ignore_dup = BooleanField(
                label="Ignore duplicate reads in read counting",
                default=False,
                description="Duplicate reads are identified using bit Ox400 in BAM FLAG "
                "field. The whole read pair is ignored if one of the reads is a "
                "duplicate read for paired-end data.",
            )

        class ExonExonJunctions:
            """Exon-exon junctions."""

            junc_counts = BooleanField(
                label="Count the number of reads supporting each exon-exon junction",
                default=False,
                description="Junctions are identified from those exon-spanning reads in "
                "input (containing 'N' in CIGAR string).",
            )

            genome = DataField(
                "seq:nucleotide",
                label="Genome",
                required=False,
                disabled="!exon_exon_junctions.junc_counts",
                description="Reference sequences used in read mapping that produced the "
                "provided BAM files. This optional argument can be used to improve read "
                "counting for junctions.",
            )

        class PairedEnd:
            """Parameters specific to paired-end reads."""

            is_paired_end = BooleanField(
                label="Count fragments (or templates) instead of reads",
                default=True,
            )

            require_both_ends_mapped = BooleanField(
                label="Count only read pairs that have both ends aligned",
                default=False,
            )

            check_frag_length = BooleanField(
                label="Check fragment length when assigning fragments to meta-features "
                "or features",
                default=False,
                description="Use minimum and maximum fragment/template length to set thresholds.",
            )

            min_frag_length = IntegerField(
                label="Minimum fragment/template length",
                default=50,
                disabled="!paired_end.check_frag_length",
            )

            max_frag_length = IntegerField(
                label="Maximum fragment/template length",
                default=600,
                disabled="!paired_end.check_frag_length",
            )

            do_not_count_chimeric_fragments = BooleanField(
                label="Do not count chimeric fragments",
                default=False,
                description="Do not count read pairs that have their two ends mapped to "
                "different chromosomes or mapped to same chromosome but on different strands.",
            )

            do_not_sort = BooleanField(
                label="Do not sort reads in BAM input",
                default=False,
            )

        class Miscellaneous:
            """Miscellaneous."""

            report_reads = BooleanField(
                label="Output detailed assignment results for each read or read pair",
                default=False,
            )

            max_mop = IntegerField(
                label="Maximum number of 'M' operations allowed in a CIGAR string",
                default=10,
                description="Both 'X' and '=' are treated as 'M' and adjacent 'M' operations "
                "are merged in the CIGAR string.",
            )

            verbose = BooleanField(
                label="Output verbose information",
                default=False,
                description="Output verbose information for debugging, such as unmatched "
                "chromosome / contig names.",
            )

        general = GroupField(General, label="General options", hidden="!advanced")

        overlap = GroupField(
            Overlap, label="Overlap between reads and features", hidden="!advanced"
        )

        read_filtering = GroupField(
            ReadFiltering, label="Read filtering", hidden="!advanced"
        )

        exon_exon_junctions = GroupField(
            ExonExonJunctions, label="Exon-exon junctions", hidden="!advanced"
        )

        paired_end = GroupField(
            PairedEnd,
            label="Parameters specific to paired-end reads",
            hidden="!advanced",
        )

        miscellaneous = GroupField(
            Miscellaneous,
            label="Miscellaneous",
            hidden="!advanced",
        )

    class Output:
        """Output fields."""

        rc = FileField(label="Read counts")
        tpm = FileField(label="TPM")
        cpm = FileField(label="CPM")
        exp = FileField(label="Normalized expression")
        exp_json = JsonField(label="Expression (json)")
        exp_type = StringField(label="Expression type")
        exp_set = FileField(label="Expressions")
        exp_set_json = JsonField(label="Expressions (json)")
        feature_counts_output = FileField(label="featureCounts output")
        counts_summary = FileField(label="Counts summary")
        read_assignments = FileField(
            label="Read assignments",
            required=False,
            description="Read assignment results for each read (or fragment if paired end).",
        )
        read_assignments = FileField(
            label="Strandedness report file",
            required=False,
        )
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        build = StringField(label="Build")
        feature_type = StringField(label="Feature type")

    def run(self, inputs, outputs):
        """Run the analysis."""

        STRANDEDNESS_CODES = {
            "IU": 0,
            "U": 0,
            "non_specific": 0,
            "ISF": 1,
            "OSF": 1,
            "SF": 1,
            "forward": 1,
            "ISR": 2,
            "OSR": 2,
            "SR": 2,
            "reverse": 2,
        }

        if inputs.aligned_reads.output.species != inputs.annotation.output.species:
            self.error(
                f"Species of aligned reads {inputs.aligned_reads.output.species} "
                f"and annotation {inputs.annotation.output.species} do not match. Please provide "
                f"aligned reads and annotation with the same species."
            )

        if inputs.aligned_reads.output.build != inputs.annotation.output.build:
            self.error(
                f"Builds of aligned reads {inputs.aligned_reads.output.species} "
                f"and annotation {inputs.annotation.output.species} do not match. Please provide "
                f"aligned reads and annotation with the same build."
            )

        if inputs.assay_type == "auto" and not inputs.cdna_index:
            self.error(
                "cDNA sequence index must be provided to automatically detect strandedness."
            )

        if (
            inputs.cdna_index
            and inputs.aligned_reads.output.species != inputs.cdna_index.output.species
        ):
            self.error(
                f"Species of aligned reads {inputs.aligned_reads.output.species} "
                f"and cDNA index {inputs.annotation.output.species} do not match. Please provide "
                f"aligned reads and cDNA index with the same species."
            )

        # Avoid reporting the full path to the alignment (.bam) file in the counts summary file.
        # This is to prevent the FeatureCounts results to be reported as a separate sample in the MultiQC report
        bam_file = Path(inputs.aligned_reads.output.bam.path).name
        Path(bam_file).symlink_to(inputs.aligned_reads.output.bam.path)

        # Set output file names
        assert bam_file.endswith(".bam")
        name = bam_file[:-4]

        exp_output = f"{name}_{inputs.normalization_type.lower()}.tab.gz"

        # check if aligned reads are single or paired-end
        paired_end = True
        if int(Cmd["samtools"]["view"]["-c", "-f", "1", bam_file]().strip()) == 0:
            paired_end = False

        self.progress = 0.05

        # set strandedness
        if inputs.assay_type == "auto":
            all_reads = int(Cmd["samtools"]["view"]["-c", bam_file]().strip())
            sampling_rate = min(inputs.n_reads / all_reads, 1)
            strand_check_bam = bam_file
            # subsample the BAM file
            if sampling_rate < 1:
                strand_check_bam = "subsampled_sorted.bam"
                (
                    Cmd["samtools"]["view"][
                        f"-@ {self.requirements.resources.cores}",
                        "-h",
                        f"-s {sampling_rate}",
                        strand_check_bam,
                    ]
                    | Cmd["samtools"]["sort"][
                        f"-@ {self.requirements.resources.cores}", "-n", "-"
                    ]
                    > strand_check_bam
                )()

            # prepare reads input
            if paired_end:
                reads_input = ["-1", "mate1.fastq", "-2", "mate2.fastq"]
            else:
                reads_input = ["-s", "reads.fastq"]

            fastq_args = [
                f"-@ {self.requirements.resources.cores}",
                "-N",
                reads_input,
                strand_check_bam,
            ]
            return_code, _, _ = Cmd["samtools"]["fastq"][fastq_args] & TEE(retcode=None)
            if return_code:
                self.error("Samtools fastq command failed.")

            salmon_out_folder = "salmon_output"

            # Run Salmon Quant
            salmon_args = [
                "-i",
                inputs.cdna_index.output.index.path,
                "-l",
                "A",
                reads_input if paired_end else ["-r", "reads.fastq"],
                "-o",
                salmon_out_folder,
                "-p",
                self.requirements.resources.cores,
                "--minAssignedFrags",
                1,
            ]
            return_code, _, _ = Cmd["salmon"]["quant"][salmon_args] & TEE(retcode=None)
            if return_code:
                self.error("Error while running Salmon Quant.")

            # Extract the strandedness code from the JSON report produced by the Salmon tool
            lib_type_report = f"{salmon_out_folder}/lib_format_counts.json"
            strand_code = json.load(open(lib_type_report)).get("expected_format", "")

            if strand_code:
                try:
                    strandedness = STRANDEDNESS_CODES[strand_code]
                except KeyError:
                    self.error(
                        f"Unsupported strand code detected: {strand_code} "
                        f"Please re-run analysis in user-selected strandedness mode."
                    )
            else:
                self.error(
                    "Automated detection of strandedness failed. "
                    "Re-run analysis in user-selected strandedness mode."
                )
        else:
            strandedness = STRANDEDNESS_CODES[inputs.assay_type]

        self.progress = 0.1

        # Replace empty gene_id entries in annotation file if source is UCSC
        annotation_file = inputs.annotation.output.annot.path
        if (
            inputs.annotation.output.source == "UCSC"
            and inputs.annotation.type.startswith("data:annotation:gtf")
        ):
            with open(annotation_file, "r") as infile:
                filedata = infile.read()

            # Replace the missing gene_ids
            annot_data = filedata.replace('gene_id "";', 'gene_id "unknown";')

            # Write the output file
            annotation_file = "annotation_modified.gtf"
            with open(annotation_file, "w") as outfile:
                outfile.write(annot_data)

        fc_output = "featureCounts_rc.txt"

        # Prepare featureCounts inputs
        args = [
            "-a",
            annotation_file,
            "-o",
            fc_output,
            "-F",
            "GTF",
            "-t",
            inputs.feature_class,
            "-g",
            inputs.id_attribute,
            "--minOverlap",
            inputs.overlap.min_overlap,
            "--fracOverlap",
            inputs.overlap.frac_overlap,
            "--fracOverlapFeature",
            inputs.overlap.frac_overlap_feature,
            "--readExtension5",
            inputs.overlap.read_extension_5,
            "--readExtension3",
            inputs.overlap.read_extension_3,
            "-Q",
            inputs.read_filtering.min_mqs,
            "--maxMOp",
            inputs.miscellaneous.max_mop,
            "-s",
            strandedness,
            "-T",
            self.requirements.resources.cores,
        ]

        if inputs.general.count_features:
            args.append("-f")

        if inputs.overlap.allow_multi_overlap:
            args.append("-O")

        if inputs.overlap.largest_overlap:
            args.append("--largestOverlap")

        if inputs.overlap.read_to_pos:
            args.extend(["--read2pos", inputs.overlap.read_to_pos])

        if inputs.general.count_multi_mapping_reads:
            args.append("-M")

        if inputs.general.fraction:
            args.append("--fraction")

        if inputs.read_filtering.split_only:
            args.append("--countSplitAlignmentsOnl")

        if inputs.read_filtering.non_split_only:
            args.append("--countNonSplitAlignmentsOnly")

        if inputs.read_filtering.primary:
            args.append("--primary")

        if inputs.read_filtering.ignore_dup:
            args.append("--ignoreDup")

        if inputs.exon_exon_junctions.junc_counts:
            args.append("-J")

        if inputs.exon_exon_junctions.junc_counts and inputs.exon_exon_junctions.genome:
            args.extend(["-G", inputs.exon_exon_junctions.genome.path])

        if inputs.general.by_read_group:
            args.append("--byReadGroup")

        if inputs.general.count_long_reads:
            args.append("-L")

        if inputs.miscellaneous.report_reads:
            args.append("-R CORE")

        if inputs.miscellaneous.verbose:
            args.append("--verbose")

        # List of options for paired-end reads
        if paired_end and inputs.paired_end.is_paired_end:
            args.append("-p")

        if paired_end and inputs.paired_end.require_both_ends_mapped:
            args.append("-B")

        if paired_end and inputs.paired_end.check_frag_length:
            args.append("-P")

        if paired_end and inputs.paired_end.check_frag_length:
            args.extend(["-d", inputs.paired_end.min_frag_length])

        if paired_end and inputs.paired_end.check_frag_length:
            args.extend(["-D", inputs.paired_end.max_frag_length])

        if paired_end and inputs.paired_end.do_not_count_chimeric_fragments:
            args.append("-C")

        if paired_end and inputs.paired_end.do_not_sort:
            args.append("--donotsort")

        # Run featureCounts
        return_code, _, _ = Cmd["featureCounts"][args][bam_file] & TEE(retcode=None)
        if return_code:
            self.error("Error while running featureCounts.")

        self.progress = 0.8

        raw_counts = "rc.txt"
        tpm = "tpm.txt"
        cpm = "cpm.txt"

        # parse featureCounts output
        get_gene_counts(infile=fc_output, outfile=raw_counts, sample_name=bam_file)
        get_gene_lenghts(infile=fc_output, outfile="gene_lengths.txt")

        # Normalize counts
        rnanorm_args = [
            raw_counts,
            "--gene-lengths",
            "gene_lengths.txt",
            "--tpm-output",
            tpm,
            "--cpm-output",
            cpm,
        ]
        return_code, _, _ = Cmd["rnanorm"][rnanorm_args] & TEE(retcode=None)
        if return_code:
            self.error("Error while normalizing counts using rnanorm.")

        self.progress = 0.9

        # prepare the expression set outputs
        feature_ids = pd.read_csv(
            raw_counts, sep="\t", index_col="FEATURE_ID"
        ).index.to_list()

        feature_filters = {
            "source": inputs.annotation.output.source,
            "species": inputs.aligned_reads.output.species,
            "feature_id__in": feature_ids,
        }

        feature_ids_to_names = {
            f.feature_id: f.name for f in self.feature.filter(**feature_filters)
        }

        prepare_expression_set(
            rc=raw_counts,
            tpm=tpm,
            cpm=cpm,
            feature_dict=feature_ids_to_names,
            outfile_name=f"{name}_expressions",
        )

        self.progress = 0.95

        # rename and compress the expression files
        rename_columns_and_compress(raw_counts, f"{name}_rc.tab.gz")
        rename_columns_and_compress(tpm, f"{name}_tpm.tab.gz")
        rename_columns_and_compress(cpm, f"{name}_cpm.tab.gz")

        # Compress the featureCounts output file
        with open(file="featureCounts_rc.txt", mode="rb") as f_in:
            with gzip.open(
                filename=f"{name}_featureCounts_rc.txt.gz", mode="wb"
            ) as f_out:
                shutil.copyfileobj(f_in, f_out)

        # Save the abundance estimates to JSON storage
        json_output = "json.txt"
        expression_to_storage(rc_input=exp_output, rc_output=json_output)

        # Save the outputs
        outputs.feature_counts_output = f"{name}_featureCounts_rc.txt.gz"
        outputs.counts_summary = "featureCounts_rc.txt.summary"
        outputs.rc = f"{name}_rc.tab.gz"
        outputs.tpm = f"{name}_tpm.tab.gz"
        outputs.cpm = f"{name}_cpm.tab.gz"
        outputs.exp = exp_output
        outputs.exp_json = json_output
        outputs.exp_set = f"{name}_expressions.txt.gz"
        outputs.exp_set_json = f"{name}_expressions.json"
        outputs.exp_type = inputs.normalization_type
        outputs.source = inputs.annotation.output.source
        outputs.species = inputs.aligned_reads.output.species
        outputs.build = inputs.aligned_reads.output.build
        outputs.feature_type = inputs.feature_type

        if inputs.miscellaneous.report_reads:
            outputs.read_assignments = f"{name}.bam.featureCounts"
