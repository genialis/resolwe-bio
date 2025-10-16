"""Import paired-end FASTQ reads."""

from pathlib import Path

from resolwe.process import (
    FileField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)

SUPPORTED_EXTENSIONS = (
    ".fastq.gz",
    ".fq.gz",
)


def rename_reads_file(mate_file):
    """Rename reads file to have .fastq.gz extension."""
    mate_file = Path(mate_file)
    extensions = "".join(mate_file.suffixes[-2:])
    renamed_reads = str(mate_file).replace(extensions, ".fastq.gz")
    mate_file.rename(renamed_reads)
    return renamed_reads


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

        # ensure that the files have the correct suffix
        for mate_file in (mate1, mate2):
            mate_file = Path(mate_file)
            if not mate_file.is_file():
                self.error("Input file {} does not exist".format(mate_file.name))
            if not mate_file.name.lower().endswith(SUPPORTED_EXTENSIONS):
                self.error(
                    "Unrecognized file name extension in file {}. "
                    "Supported file name extensions are {}.".format(
                        mate_file.name, SUPPORTED_EXTENSIONS
                    )
                )

        mate1 = rename_reads_file(mate1)
        mate2 = rename_reads_file(mate2)

        # save the outputs
        outputs.mate1 = mate1
        outputs.mate2 = mate2

        # set sample species annotation
        if inputs.species:
            outputs.species = inputs.species
            self.data.entity.annotations["general.species"] = inputs.species
