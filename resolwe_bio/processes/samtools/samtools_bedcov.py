"""Samtools bedcov."""

from io import StringIO
from pathlib import Path

import pandas as pd
from plumbum import TEE

from resolwe.process import (
    BooleanField,
    Cmd,
    DataField,
    FileField,
    GroupField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)


def normalize_bedcov_output(output: str) -> pd.DataFrame:
    """Calculate mean coverage over regions specified in samtools bedcov output.

    Samtools bedcov returns sum of read coverages across all bases in bed file regions.
    This function normalizes the values with the length of individual regions.

    Input columns : [chr, start, stop, coverage sum]
    Output columns : [chr, start, stop, coverage mean]
    """
    coverage_df = pd.read_csv(
        StringIO(output),
        sep="\t",
        names=["chrom", "start", "end", "coverage"],
    )

    coverage_df["coverage"] = coverage_df["coverage"] / (
        coverage_df["end"] - coverage_df["start"]
    )
    coverage_df["coverage"] = coverage_df["coverage"].round(2)

    return coverage_df


class SamtoolsBedcov(Process):
    """Samtools bedcov.

    Reports the total read base count (i.e. the sum of per base read depths)
    for each genomic region specified in the supplied BED file.
    The regions are output as they appear in the BED file and are 0-based.
    The output is tab-delimited with the first three columns being the chromosome,
    start and end positions of the region, and the remaining column
    being the sum read base counts or normalized sum read base counts by the length of
    individual regions (mean coverage).

    For more information about samtools bedcov, click
    [here](https://www.htslib.org/doc/samtools-bedcov.html).

    """

    slug = "samtools-bedcov"
    process_type = "data:alignment:bam:samtools"
    name = "Samtools bedcov"
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
    version = "1.0.1"
    entity = {"type": "sample"}
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields for SamtoolsBedcov."""

        bam = DataField(data_type="alignment:bam", label="Input BAM file")
        bedfile = DataField(
            data_type="bed",
            label="Target BED file",
            description="Target BED file with regions to extract.",
        )

        class AdvancedOptions:
            """Advanced options."""

            min_read_qual = IntegerField(
                label="Only count reads with mapping quality greater than or equal to [-Q]",
                required=False,
            )
            rm_del_ref_skips = BooleanField(
                label="Do not include deletions (D) and ref skips (N) in bedcov computation. [-j]",
                default=False,
            )
            output_option = StringField(
                label="Metric by which to output coverage.",
                default="sum",
                choices=[
                    ("sum", "Sum (default)"),
                    ("mean", "Mean"),
                ],
                required=False,
                description="Choose to output the sum of reads per base of specified regions,"
                "or these values normalized by their length."
                "This method is not part of the bedcov command, but is a custom option added to the process.",
            )

        advanced = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for SamtoolsBedcov."""

        coverage_report = FileField(label="Output coverage report")

    def run(self, inputs, outputs):
        """Run the analysis."""

        name = f"{Path(inputs.bedfile.output.bed.path).stem}"
        coverage_report = f"{name}_coverage.txt"

        if inputs.bedfile.output.species != inputs.bam.output.species:
            self.error(
                "Input BAM file and BED file are of different species. "
                f"BAM file is from {inputs.bam.output.species}, "
                f"while BED file is from {inputs.bedfile.output.species}."
            )
        if inputs.bedfile.output.build != inputs.bam.output.build:
            self.error(
                "Input BAM file and BED file have different genome build, "
                "but it should be the same. BAM file has build "
                f"{inputs.bam.output.build}, while BED file has build "
                f"{inputs.bedfile.output.build}."
            )

        input_options = []

        if inputs.advanced.min_read_qual:
            input_options.extend(["-Q", inputs.advanced.min_read_qual])
        if inputs.advanced.rm_del_ref_skips:
            input_options.append("-j")

        cmd = Cmd["samtools"]["bedcov"][input_options][
            inputs.bedfile.output.bed.path,
            inputs.bam.output.bam.path,
        ]
        return_code, stdout, stderr = cmd & TEE(retcode=None)
        if return_code:
            self.error(f"Samtools bedcov failed. {stdout}, {stderr}")

        if inputs.advanced.output_option == "mean":
            coverage_df = normalize_bedcov_output(stdout)
            coverage_df.to_csv(
                coverage_report,
                sep="\t",
                index=False,
                header=False,
            )
        else:
            with open(coverage_report, "w") as f:
                f.writelines(stdout)

        outputs.coverage_report = coverage_report
