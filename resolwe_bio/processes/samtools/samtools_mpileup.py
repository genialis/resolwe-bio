"""Samtools mpileup."""

import re
from pathlib import Path

import numpy as np
import pandas as pd
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    ListField,
    Process,
    SchedulingClass,
    StringField,
)

INDEL_PATTERN = r"([+-])(\d+)(\w+)"


def remove_indel_chars(match):
    """Remove adjacent indel base chars in the sequence, leaving just + or -."""
    indel_length = int(match.group(2))
    indel_string = match.group(3)
    return match.group(1) + indel_string[indel_length:]


def prepare_char_counts(pileup_output, column_names, self, remove_read_ends):
    """Generate a dataframe with character counts."""
    df = pd.read_csv(
        pileup_output,
        sep="\t",
        names=column_names,
    )
    if df.empty:
        self.error("No data in the pileup file.")

    pos_df = df[["chr", "pos"]].copy()

    char_count_columns = [
        "A",
        "C",
        "T",
        "G",
        "N",
        "-",
        "+",
        "<",
        ">",
        "*",
    ]
    if not remove_read_ends:
        char_count_columns.extend(["^", "$"])

    char_counts = np.zeros((len(df), len(char_count_columns)), dtype=int)

    ref = df["ref"].values
    sequences = df["base_sequence"].values

    for i, sequence in enumerate(sequences):
        sequence = str(sequence).upper()

        # Remove mapping quality after the ^ character
        if not remove_read_ends:
            sequence = re.sub(r"\^(.)", "^", sequence)

        # Delete indel sequences, leaving just +/-
        sequence = re.sub(INDEL_PATTERN, remove_indel_chars, sequence)

        ref_base = ref[i].upper()

        # Handle counts for reference matches
        char_counts[i, char_count_columns.index(ref_base)] += sequence.count(
            "."
        ) + sequence.count(",")

        # Handle counts for bases and special characters
        for j, char in enumerate(char_count_columns):
            char_counts[i, j] += sequence.count(char)

    base_counts_df = pd.DataFrame(char_counts, columns=char_count_columns)

    out_df = pd.concat([pos_df, base_counts_df], axis=1)
    return out_df


class SamtoolsMpileupSingle(Process):
    """Samtools mpileup for a single BAM file.

    For more information about samtools mpileup, click
    [here](https://www.htslib.org/doc/samtools-mpileup.html).

    """

    slug = "samtools-mpileup-single"
    process_type = "data:samtoolsmpileup:single"
    name = "Samtools mpileup (single-sample)"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"},
        },
        "resources": {
            "cores": 2,
            "memory": 8192,
        },
    }
    category = "Samtools"
    data_name = "{{ bam|name|default('?') }}"
    version = "2.0.0"
    entity = {"type": "sample"}
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields for SamtoolsMpileupSingle."""

        bam = DataField(data_type="alignment:bam", label="Input BAM file")
        bedfile = DataField(
            data_type="bed",
            label="BED file",
            description="BED file containing a list of regions or sites where pileup should be generated. "
            "BED files contain at least 3 columns (chromosome, start and end position) and are 0-based half-open. [-l]",
            required=False,
            disabled="positions",
        )
        positions = DataField(
            data_type="file",
            label="List of positions",
            description="File containing a list of regions or sites where pileup should be generated. "
            "Position list files contain two columns (chromosome and position) and start counting from 1. [-l]",
            required=False,
            disabled="bedfile",
        )
        fasta = DataField(
            data_type="seq:nucleotide",
            label="Reference genome sequence",
            description="Reference genome sequence in FASTA format. "
            "It is used to display the reference base in the pileup output and do base alignment quality (BAQ) computation. [-f]",
            required=False,
        )

        class AdvancedOptions:
            """Advanced options."""

            disable_baq = BooleanField(
                label="Disable BAQ computation",
                description="BAQ is turned on when a reference file is supplied. "
                "This option disables BAQ computation. [-B]",
                default=False,
            )
            ignore_overlap_removal = BooleanField(
                label="Ignore overlap removal",
                description="Ignore the removal of overlaps."
                "When enabled, where the ends of a read-pair overlap the overlapping region will have one base selected "
                "and the duplicate base nullified by setting its phred score to zero. "
                "It will then be discarded by the --min-BQ option unless this is zero. [-x]",
                default=False,
            )
            region = StringField(
                label="Region",
                description="Region to include in the pileup. "
                "If used in conjunction with a BED file or List of positions, it then considers the intersection of the two requests. "
                "Region should be specified in the format 'chr:start-end'.",
                required=False,
            )
            max_depth = IntegerField(
                label="Max depth",
                description="Maximum depth to avoid excessive memory usage. "
                "Setting this limit reduces the amount of memory and time needed to process regions with very high coverage. "
                "Passing zero for this option sets it to the highest possible value, effectively removing the depth limit. [-d]",
                default=0,
            )
            output_mapping_qual = BooleanField(
                label="Output mapping quality",
                description="Output mapping quality as a separate column in the pileup output. [-s]",
                default=False,
            )
            output_read_names = BooleanField(
                label="Output read names",
                description="Output read names in the pileup output. [--output-QNAME]",
                default=False,
            )
            output_zero_depth = BooleanField(
                label="Output all positions (including those with zero depth)",
                description="Output all positions including those with zero depth. [-a]",
                default=False,
                disabled="advanced.output_all_positions == true",
            )
            output_all_positions = BooleanField(
                label="Output absolutely all positions",
                description="Output absolutely all positions, including unused reference sequences. "
                "Note that when used in conjunction with a BED file the Output all positions option may "
                "sometimes operate as if Output absolutely all positions was specified. This can happen "
                "if the reference sequence has coverage outside of the region specified in the BED file. [-aa]",
                default=False,
            )
            min_base_qual = IntegerField(
                label="Minimum base quality for a base to be considered",
                description="Only count reads with base quality greater than or equal to value specific. [-Q]",
                required=False,
            )
            min_mapping_qual = IntegerField(
                label="Minimum mapping quality for an alignment to be used",
                description="Only count reads with mapping quality greater than or equal to value specific. [-q]",
                required=False,
            )
            excl_flags = ListField(
                StringField(),
                label="Filter flags",
                default=["SECONDARY", "QCFAIL", "DUP"],
                description="Filter flags: skip reads with mask bits set. "
                "Press ENTER after each flag. [-ff]",
            )
            excl_output_ends = BooleanField(
                label="Exclude output ends in sequence column",
                description="Removes the “^” (with mapping quality) and “$” markup from the sequence column. [--no-output-ends]",
                default=False,
            )
            output_base_counts = BooleanField(
                label="Output per base counts",
                description="Output base counts at each position. "
                "This is a custom implementation and not part of the original samtools mpileup. "
                "It outputs a file with base counts at each position, including ref skips, delins and read ends.",
                default=True,
            )

        advanced = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for SamtoolsMpileupSingle."""

        pileup_report = FileField(label="Output pileup report")
        base_counts = FileField(label="Output per base counts", required=False)
        species = StringField(label="Species")
        build = StringField(label="Genome build")

    def run(self, inputs, outputs):
        """Run the analysis."""

        name = Path(inputs.bam.output.bam.path).stem
        result_fn = f"{name}_mpileup.txt"

        if inputs.positions and inputs.bedfile:
            self.error(
                "Only one of BED file or List of positions can be provided as input."
            )

        input_options = []
        column_names = [
            "chr",
            "pos",
            "ref",
            "depth",
            "base_sequence",
            "base_qualities",
        ]

        if inputs.fasta:
            if inputs.fasta.output.species != inputs.bam.output.species:
                self.error(
                    "Input BAM file and FASTA file are of different species. "
                    f"BAM file is from {inputs.bam.output.species}, "
                    f"while FASTA file is from {inputs.fasta.output.species}."
                )
            if inputs.fasta.output.build != inputs.bam.output.build:
                self.error(
                    "Input BAM and FASTA files should have matching genome builds. "
                    "BAM file has build "
                    f"{inputs.bam.output.build}, while FASTA file has build "
                    f"{inputs.fasta.output.build}."
                )
            input_options.extend(["-f", inputs.fasta.output.fasta.path])

        if inputs.bedfile:
            if inputs.bedfile.output.species != inputs.bam.output.species:
                self.error(
                    "Input BAM file and BED file are of different species. "
                    f"BAM file is from {inputs.bam.output.species}, "
                    f"while BED file is from {inputs.bedfile.output.species}."
                )
            if inputs.bedfile.output.build != inputs.bam.output.build:
                self.error(
                    "Input BAM and BED files should have matching genome builds. "
                    "BAM file has build "
                    f"{inputs.bam.output.build}, while BED file has build "
                    f"{inputs.bedfile.output.build}."
                )
            input_options.extend(["-l", inputs.bedfile.output.bed.path])
        elif inputs.positions:
            try:
                positions_df = pd.read_csv(
                    inputs.positions.output.file.path, sep="\t", header=None
                )
                if not positions_df.shape[1] == 2:
                    self.error(
                        "The file does not appear to be a valid list of positions. "
                        "It should contain two columns (chromosome and position) and be tab separated."
                    )
            except pd.errors.ParserError:
                self.error("Check positions file integrity.")
            input_options.extend(["-l", inputs.positions.output.file.path])

        input_options.extend(["-d", inputs.advanced.max_depth])
        if inputs.advanced.output_zero_depth:
            input_options.append("-a")
        if inputs.advanced.output_all_positions:
            input_options.append("-aa")
        if inputs.advanced.disable_baq:
            input_options.append("-B")
        if inputs.advanced.ignore_overlap_removal:
            input_options.append("-x")
        if inputs.advanced.min_base_qual:
            input_options.extend(["-Q", inputs.advanced.min_base_qual])
        if inputs.advanced.min_mapping_qual:
            input_options.extend(["-q", inputs.advanced.min_mapping_qual])
        if inputs.advanced.excl_output_ends:
            input_options.append("--no-output-ends")
        if inputs.advanced.excl_flags:
            flags = ",".join(inputs.advanced.excl_flags)
            input_options.extend(["--ff", flags])
        if inputs.advanced.region:
            input_options.extend(["-r", inputs.advanced.region])
        if inputs.advanced.output_mapping_qual:
            input_options.append("-s")
            column_names.append("mapping_quality")
        if inputs.advanced.output_read_names:
            input_options.append("--output-QNAME")
            column_names.append("read_names")

        return_code, stdout, _ = Cmd["samtools"]["mpileup"][input_options][
            inputs.bam.output.bam.path
        ]["-o", result_fn] & TEE(retcode=None)
        if return_code:
            self.error("Samtools mpileup failed.")
        if "[E::faidx_adjust_position]" in stdout:
            self.warning(
                "Samtools mpileup returned a warning due to an issue with the reference genome. "
                "Please check that the reference genome is correctly formatted and indexed."
            )

        if inputs.advanced.output_base_counts:
            depth_df = prepare_char_counts(
                pileup_output=result_fn,
                column_names=column_names,
                self=self,
                remove_read_ends=inputs.advanced.excl_output_ends,
            )
            base_counts_fn = f"{name}_base_counts.txt.gz"
            depth_df.to_csv(base_counts_fn, sep="\t", index=False, compression="gzip")
            outputs.base_counts = base_counts_fn

        result_fn_zipped = result_fn + ".gz"

        return_code, _, _ = Cmd["gzip"][result_fn] & TEE(retcode=None)
        if return_code:
            self.error("Compression of file failed.")

        outputs.pileup_report = result_fn_zipped
        outputs.species = inputs.bam.output.species
        outputs.build = inputs.bam.output.build
