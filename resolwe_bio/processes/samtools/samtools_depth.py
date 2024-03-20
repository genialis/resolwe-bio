"""Samtools depth."""

import gzip
import os
from pathlib import Path

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


class SamtoolsDepthSingle(Process):
    """Samtools depth for single BAM file.

    Computes the depth at each position or region.

    For more information about samtools depth, click
    [here](https://www.htslib.org/doc/samtools-depth.html).

    """

    slug = "samtools-depth-single"
    process_type = "data:samtoolsdepth:single"
    name = "Samtools depth (single-sample)"
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
    version = "1.1.0"
    entity = {"type": "sample"}
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields for SamtoolsDepthSingle."""

        bam = DataField(data_type="alignment:bam", label="Input BAM file")
        bedfile = DataField(
            data_type="bed",
            label="BED file",
            description="Target BED file with regions to extract. If not selected, the whole genome will be used.",
            required=False,
        )

        class AdvancedOptions:
            """Advanced options."""

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
            ignore_reads = IntegerField(
                label="Ignore reads shorter than specified length",
                description="Input the specified length under which reads are ignored. "
                "This is the number of bases in the sequence, minus any soft clips. [-I]",
                required=False,
            )
            min_base_qual = IntegerField(
                label="Minimum base quality",
                description="Only count reads with base quality greater than or equal to value specific. [-q]",
                required=False,
            )
            min_mapping_qual = IntegerField(
                label="Minimum read mapping quality",
                description="Only count reads with mapping quality greater than or equal to value specific. [-Q]",
                required=False,
            )
            excl_flags = ListField(
                StringField(),
                label="Filter flags",
                default=["UNMAP", "SECONDARY", "QCFAIL", "DUP"],
                description="Filter flags: skip reads with mask bits set. "
                "Press ENTER after each flag. [-G]",
            )
            include_deletions = BooleanField(
                label="Include deletions",
                description="Include reads with deletions in the depth computation. [-J]",
                default=False,
            )
            collapse_overlaps = BooleanField(
                label="Collapse overlapping read pairs",
                description="For the overlapping section of a read pair, count only the bases of the first read. [-s]",
                default=False,
            )
            prune_zero_depth = BooleanField(
                label="Prune zero depth entries",
                description="Prune zero depth entries from the output.",
                default=False,
            )

        advanced = GroupField(AdvancedOptions, label="Advanced options")

    class Output:
        """Output fields for SamtoolsDepthSingle."""

        depth_report = FileField(label="Output depth report")
        species = StringField(label="Species")
        build = StringField(label="Genome build")

    def run(self, inputs, outputs):
        """Run the analysis."""

        name = Path(inputs.bam.output.bam.path).stem
        result_fn = f"{name}_depth.txt"

        input_options = []

        if inputs.advanced.output_zero_depth:
            input_options.append("-a")
        if inputs.advanced.output_all_positions:
            input_options.append("-aa")
        if inputs.advanced.ignore_reads:
            input_options.extend(["-I", inputs.advanced.ignore_reads])
        if inputs.advanced.min_base_qual:
            input_options.extend(["-q", inputs.advanced.min_base_qual])
        if inputs.advanced.min_mapping_qual:
            input_options.extend(["-Q", inputs.advanced.min_mapping_qual])
        if inputs.advanced.excl_flags:
            flags = ",".join(inputs.advanced.excl_flags)
            input_options.extend(["-G", flags])
        if inputs.advanced.include_deletions:
            input_options.append("-J")
        if inputs.advanced.collapse_overlaps:
            input_options.append("-s")
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
            input_options.extend(["-b", inputs.bedfile.output.bed.path])

        return_code, stdout, stderr = Cmd["samtools"]["depth"][input_options][
            inputs.bam.output.bam.path
        ]["-o", result_fn] & TEE(retcode=None)
        if return_code:
            self.error(f"Samtools depth failed. {stdout}, {stderr}")

        result_fn_zipped = result_fn + ".gz"

        if inputs.advanced.prune_zero_depth:
            with open(result_fn, "r") as infile, gzip.open(
                result_fn_zipped, "w"
            ) as outfile:
                for line in infile:
                    columns = line.strip().split("\t")
                    if columns[2] != "0":
                        outfile.write(line.encode())

            os.remove(result_fn)

        else:
            return_code, _, _ = Cmd["gzip"][result_fn] & TEE(retcode=None)
            if return_code:
                self.error("Compression of file failed.")

        outputs.depth_report = result_fn_zipped
        outputs.species = inputs.bam.output.species
        outputs.build = inputs.bam.output.build
