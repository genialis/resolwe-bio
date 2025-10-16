"""Import paired-end FASTQ reads."""
from pathlib import Path

from resolwe.process import (
    FileField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


def ensure_fastq_gz(file_path: str) -> str:
    """Ensure the file has a .fastq.gz suffix by renaming if needed."""
    path = Path(file_path)
    suffixes = [suffix.lower() for suffix in path.suffixes]
    if suffixes[-2:] == [".fastq", ".gz"]:
        return str(path)

    base = path
    for _ in path.suffixes:
        base = base.with_suffix("")

    new_path = base.with_name(f"{base.name}.fastq.gz")

    if path != new_path:
        path.rename(new_path)

    return str(new_path)


class UploadFastqPaired(Process):
    """Import paired-end reads in FASTQ format."""

    slug = "upload-fastq-paired-docs"
    name = "FASTQ file (paired-end)"
    process_type = "data:reads:fastq:paired"
    version = "1.0.0"
    category = "Import"
    data_name = '{{ mate1.file|default("?") }}'
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    entity = {
        "type": "sample",
    }
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"}
        },
        "resources": {
            "cores": 1,
            "memory": 2048,
            "storage": 100,
            "network": False,
        },
    }

    class Input:
        """Input fields to process UploadFastqPaired."""

        mate1 = FileField(
            label="Mate 1",
            description="Sequencing reads in FASTQ format. Mate 1 file."
        )
        mate2 = FileField(
            label="Mate 2",
            description="Sequencing reads in FASTQ format. Mate 2 file."
        )
        species = StringField(
            label="Species",
            description="Species of the sample.",
            required=False,
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
            ]
        )


    class Output:
        """Output fields to process UploadFastqPaired."""

        mate1 = FileField(label="Reads file (mate 1)")
        mate2 = FileField(label="Reads file (mate 2)")
        species = StringField(label="Species", required=False)

    def run(self, inputs, outputs):
        """Run upload."""

        # import FASTQ files in compressed format
        mate1 = inputs.mate1.import_file(imported_format="compressed")
        mate2 = inputs.mate2.import_file(imported_format="compressed")

        # ensure that the files have the correct suffix and save to the
        # output fields
        outputs.mate1 = ensure_fastq_gz(mate1)
        outputs.mate2 = ensure_fastq_gz(mate2)

        # set sample species annotation
        if inputs.species:
            outputs.species = inputs.species
            self.data.entity.annotations["general.species"] = inputs.species
