"""Converting BAM to BEDPE and normalized BigWig files."""

import os
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    FloatField,
    IntegerField,
    Process,
    SchedulingClass,
    StringField,
)


class BamToBedpe(Process):
    """Takes in a BAM file and calculates a normalization factor in BEDPE format.

    Done by sorting with Samtools and transformed with Bedtools.
    """

    slug = "bedtools-bamtobed"
    name = "Bedtools bamtobed"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {"cores": 1, "memory": 8192},
    }
    data_name = "{{ alignment|name|default('?') }}"
    version = "1.3.1"
    process_type = "data:bedpe"
    category = "BAM processing"
    entity = {"type": "sample"}
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        alignment = DataField("alignment:bam", label="Alignment BAM file")

    class Output:
        """Output fields."""

        bedpe = FileField(label="BEDPE file")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""
        path = inputs.alignment.output.bam.path
        basename = os.path.basename(path)
        assert basename.endswith(".bam")
        name = basename[:-4]
        bedpe_file = f"{name}.bedpe"

        samtools_param = ["-n", path]
        bedtools_param = ["-bedpe", "-i"]

        (
            Cmd["samtools"]["sort"][samtools_param]
            | Cmd["bedtools"]["bamtobed"][bedtools_param]
            > bedpe_file
        )()

        if not os.path.exists(bedpe_file):
            self.error("Converting BAM to BEDPE with Bedtools bamtobed failed.")

        outputs.bedpe = bedpe_file
        outputs.species = inputs.alignment.output.species
        outputs.build = inputs.alignment.output.build


class CalculateBigWig(Process):
    """Calculate bigWig coverage track.

    Deeptools bamCoverage takes an alignment of reads or fragments as
    input (BAM file) and generates a coverage track (bigWig) as output.
    The coverage is calculated as the number of reads per bin, where
    bins are short consecutive counting windows of a defined size. For
    more information is available in the
    [bamCoverage documentation](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html).
    """

    slug = "calculate-bigwig"
    name = "Calculate coverage (bamCoverage)"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {"cores": 1, "memory": 16384},
    }
    data_name = "{{ alignment|name|default('?') }}"
    version = "2.0.1"
    process_type = "data:coverage:bigwig"
    category = "BAM processing"
    entity = {"type": "sample"}
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input fields."""

        alignment = DataField("alignment:bam", label="Alignment BAM file")
        bedpe = DataField(
            "bedpe",
            label="BEDPE Normalization factor",
            description="The BEDPE file describes disjoint genome features, "
            "such as structural variations or paired-end sequence alignments. "
            "It is used to estimate the scale factor [--scaleFactor].",
            required=False,
        )
        scale = FloatField(
            label="Scale for the normalization factor",
            description="Magnitude of the scale factor. The scaling factor "
            "[--scaleFactor] is calculated by dividing the scale with the "
            "number of features in BEDPE (scale/(number of features)).",
            disabled="!bedpe",
            default=10000,
        )
        bin_size = IntegerField(
            label="Bin size[--binSize]",
            description="Size of the bins (in bp) for the output bigWig file. "
            "A smaller bin size value will result in a higher resolution of "
            "the coverage track but also in a larger file size.",
            default=50,
        )

    class Output:
        """Output fields."""

        bigwig = FileField(label="Coverage file (bigWig)")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""
        bam_path = Path(inputs.alignment.output.bam.path)
        assert bam_path.name.endswith(".bam")
        name = bam_path.stem

        if inputs.bedpe:
            with open(inputs.bedpe.output.bedpe.path, "rb") as f:
                spike_count = sum(1 for _ in f)

            if spike_count == 0:
                self.error("BEDPE file is empty there were no features found.")

            scale_factor = inputs.scale / spike_count
            out_file = Path(f"{name}.SInorm.bigwig")
        else:
            scale_factor = 1
            out_file = Path(f"{name}.bigwig")

        self.progress(0.1)

        bam_coverage_param = [
            "--bam",
            bam_path,
            "--scaleFactor",
            scale_factor,
            "--outFileName",
            out_file,
            "--numberOfProcessors",
            self.requirements.resources.cores,
            "--outFileFormat",
            "bigwig",
            "--binSize",
            inputs.bin_size,
        ]

        return_code, _, stderr = Cmd["bamCoverage"][bam_coverage_param] & TEE(
            retcode=None
        )

        if return_code:
            print(stderr)
            self.error("Calculating coverage with bamCoverage failed.")

        if not out_file.is_file():
            self.error("Generation of a scaled bigWig file failed.")

        self.progress(0.9)

        outputs.bigwig = str(out_file)
        outputs.species = inputs.alignment.output.species
        outputs.build = inputs.alignment.output.build
