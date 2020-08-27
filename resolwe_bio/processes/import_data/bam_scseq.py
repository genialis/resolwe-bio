"""Upload single cell BAM."""
import os
from shutil import move

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
    version = "1.1.1"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    entity = {"type": "sample"}
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/common:1.3.1"}},
    }
    data_name = '{{ reads|sample_name|default("?") }}'

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
        bam_path = os.path.basename(inputs.src.path)
        bai_path = os.path.basename(inputs.src2.path)
        assert bam_path.endswith(".bam")
        assert bai_path.endswith(".bam.bai")
        bam_name = bam_path[:-4]
        bai_name = bai_path[:-8]
        if bam_name != bai_name:
            self.error("BAM and BAI files should have the same name.")

        move(inputs.src.file_temp, bam_path)
        move(inputs.src2.file_temp, bai_path)

        stats = "{}_stats.txt".format(bam_name)
        (Cmd["samtools"]["flagstat"][bam_path] > stats)()

        outputs.bam = bam_path
        outputs.bai = bai_path
        outputs.stats = stats
        outputs.species = inputs.species
        outputs.build = inputs.build
