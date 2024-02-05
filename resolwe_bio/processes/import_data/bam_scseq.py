"""Upload single cell BAM."""

from resolwe.process import (
    Cmd,
    DataField,
    FileField,
    Process,
    SchedulingClass,
    StringField,
)


class ImportScBam(Process):
    """Import scSeq BAM file and index."""

    slug = "upload-bam-scseq-indexed"
    name = "Single cell BAM file and index"
    process_type = "data:alignment:bam:scseq"
    version = "1.4.1"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/common:4.1.1"}
        },
    }
    data_name = "{{ reads|name|default('?') }}"

    class Input:
        """Input fields to process Import ScBam."""

        src = FileField(
            description="A mapping file in BAM format.",
            label="Mapping (BAM)",
        )
        src2 = FileField(
            description="An index file of a BAM mapping file (ending with bam.bai).",
            label="BAM index (*.bam.bai file)",
        )
        reads = DataField(
            data_type="screads:",
            label="Single cell fastq reads",
        )
        species = StringField(
            label="Species",
            description="Species latin name.",
        )
        build = StringField(
            label="Build",
        )

    class Output:
        """Output fields to process Import ScBam."""

        bam = FileField(label="Uploaded BAM")
        bai = FileField(label="Index BAI")
        stats = FileField(label="Alignment statistics")
        build = StringField(label="Build")
        species = StringField(label="Species")

    def run(self, inputs, outputs):
        """Run the analysis."""
        bam_path = inputs.src.import_file(imported_format="extracted")
        bai_path = inputs.src2.import_file(imported_format="extracted")
        assert bam_path.endswith(".bam")
        assert bai_path.endswith(".bam.bai")
        bam_name = bam_path[:-4]
        bai_name = bai_path[:-8]
        if bam_name != bai_name:
            self.error("BAM and BAI files should have the same name.")

        stats = "{}_stats.txt".format(bam_name)
        (Cmd["samtools"]["flagstat"][bam_path] > stats)()

        outputs.bam = bam_path
        outputs.bai = bai_path
        outputs.stats = stats
        outputs.species = inputs.species
        outputs.build = inputs.build
