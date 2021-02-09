"""Converting BAM to BEDPE and normalized BigWig files."""
import os

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    FloatField,
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
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/rnaseq:5.12.0"}
        },
        "resources": {"cores": 1, "memory": 8192},
    }
    data_name = "Bedtools bamtobed ({{alignment|sample_name|default('?')}})"
    version = "1.1.1"
    process_type = "data:bedpe"
    category = "Other"
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


class ScaleBigWig(Process):
    """Creates a scaled BigWig file."""

    slug = "scale-bigwig"
    name = "Deeptools bamCoverage"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/rnaseq:5.12.0"}
        },
        "resources": {"cores": 1, "memory": 16384},
    }
    data_name = "Scale BigWig ({{alignment|sample_name|default('?')}})"
    version = "1.1.1"
    process_type = "data:coverage:bigwig"
    category = "Other"
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
            "It is used to estimate the scale factor.",
        )
        scale = FloatField(
            label="Scale for the normalization factor",
            description="Magnitude of the scale factor. "
            "The scaling factor is calculated by dividing the scale "
            "with the number of features in BEDPE "
            "(scale/(number of features)).",
            default=10000,
        )

    class Output:
        """Output fields."""

        bigwig = FileField(label="bigwig file")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""
        path = inputs.alignment.output.bam.path
        basename = os.path.basename(path)
        assert basename.endswith(".bam")
        name = basename[:-4]
        out_file = f"{name}.SInorm.bigwig"
        out_index = f"{name}.bai"

        with open(inputs.bedpe.output.bedpe.path) as f:
            spike_count = f.readlines()
            spike_count = len(spike_count)
        scale_factor = inputs.scale / spike_count

        bam_coverage_param = [
            "--bam",
            path,
            "--scaleFactor",
            scale_factor,
            "--outFileName",
            out_file,
            "--numberOfProcessors",
            self.requirements.resources.cores,
            "--outFileFormat",
            "bigwig",
        ]

        (Cmd["samtools"]["index"][path][out_index])()
        self.progress(0.5)
        (Cmd["bamCoverage"][bam_coverage_param])()

        if not os.path.exists(out_file):
            self.error("Generation of a scaled BigWig file with bamCoverage failed.")

        outputs.bigwig = out_file
        outputs.species = inputs.alignment.output.species
        outputs.build = inputs.alignment.output.build
