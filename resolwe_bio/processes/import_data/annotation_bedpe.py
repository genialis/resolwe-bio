"""Upload BEDPE files."""

from resolwe.process import FileField, Process, SchedulingClass, StringField


class ImportBEDPEFile(Process):
    """Upload BEDPE files."""

    slug = "upload-bedpe"
    name = "BEDPE file"
    process_type = "data:bedpe:"
    data_name = '{{ src.file|default("?") }}'
    version = "1.1.0"
    category = "Import"
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/common:1.3.1",},},
        "resources": {"cores": 1, "memory": 1024, "network": True,},
    }
    scheduling_class = SchedulingClass.BATCH

    class Input:
        """Input parameters."""

        src = FileField(label="Select BEDPE file to upload")
        species = StringField(
            label="Species",
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
                ("Odocoileus virginianus texanus", "Odocoileus virginianus texanus"),
                ("Solanum tuberosum", "Solanum tuberosum"),
            ],
        )
        build = StringField(label="Build")

    class Output:
        """Output parameters."""

        bedpe = FileField(label="BEDPE file")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Upload BEDPE file."""
        # This process consists of upload and zipping process.
        # A validator of this format would be nice.
        # https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
        bedpe_path = inputs.src.import_file()

        outputs.bedpe = bedpe_path
        outputs.species = inputs.species
        outputs.build = inputs.build
