"""Upload geneset."""

import gzip
import json
from pathlib import Path

from resolwe.process import (
    FileField,
    JsonField,
    ListField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)
from resolwe.process.models import DescriptorSchema


def parse_geneset_file(geneset_file, warning):
    """Parse geneset file."""
    with open(geneset_file, "rU") as handle:
        # skip empty lines
        genes = [str(line.strip()) for line in handle if line.strip()]
        geneset = sorted(set(genes))
        if len(genes) != len(geneset):
            warning("Removed duplicated genes.")

    return geneset


def save_geneset_to_json(geneset, output_json):
    """Save geneset to json file."""
    with open(output_json, "w") as handle:
        json.dump(
            {"genes": geneset},
            handle,
            separators=(",", ":"),
            allow_nan=False,
        )


def save_geneset_to_file(geneset, output_file):
    """Save geneset to file."""
    with gzip.open(output_file, "w") as handle:
        handle.write("\n".join(geneset).encode("utf-8"))


class UploadGeneset(Process):
    """
    Upload a set of genes.

    Provide one gene ID per line in a .tab, .tab.gz, or .txt file format.
    """

    slug = "upload-geneset"
    name = "Gene set"
    process_type = "data:geneset"
    version = "1.3.2"
    category = "Import"
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.RAW
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/genialis/resolwebio/base:ubuntu-22.04-14112023"
            }
        },
        "resources": {
            "cores": 1,
            "memory": 1024,
            "storage": 10,
        },
    }
    data_name = '{{ src.file|default("?") }}'

    class Input:
        """Input fields."""

        src = FileField(
            label="Gene set",
            description="List of genes (.tab/.txt extension), one gene ID per line.",
        )
        source = StringField(
            label="Gene ID source",
            allow_custom_choice=True,
            choices=[
                ("AFFY", "AFFY"),
                ("DICTYBASE", "DICTYBASE"),
                ("ENSEMBL", "ENSEMBL"),
                ("NCBI", "NCBI"),
                ("UCSC", "UCSC"),
            ],
        )
        species = StringField(
            label="Species",
            description="Species latin name.",
            allow_custom_choice=True,
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
                ("Odocoileus virginianus texanus", "Odocoileus virginianus texanus"),
                ("Solanum tuberosum", "Solanum tuberosum"),
            ],
        )

    class Output:
        """Output fields."""

        geneset = FileField(label="Gene set")
        geneset_json = JsonField(label="Gene set (JSON)")
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")

    def run(self, inputs, outputs):
        """Run the analysis."""
        geneset_file = inputs.src.import_file(imported_format="extracted")
        geneset = parse_geneset_file(geneset_file, self.warning)

        # Save geneset to file
        save_geneset_to_file(geneset, f"{Path(geneset_file).stem}.tab.gz")
        outputs.geneset = f"{Path(geneset_file).stem}.tab.gz"

        # Save geneset to json
        save_geneset_to_json(geneset, "geneset.json")
        outputs.geneset_json = "geneset.json"

        outputs.source = inputs.source
        outputs.species = inputs.species

        # Set the descriptor schema of this object to geneset
        ds = DescriptorSchema.get_latest(slug="geneset")
        self.data.descriptor_schema = ds.id


class CreateGeneset(Process):
    """Create a gene set from a list of genes."""

    slug = "create-geneset"
    name = "Gene set (create)"
    process_type = "data:geneset"
    version = "1.3.2"
    category = "Import"
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.RAW
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/genialis/resolwebio/base:ubuntu-22.04-14112023"
            }
        },
        "resources": {
            "cores": 1,
            "memory": 1024,
            "storage": 10,
        },
    }
    data_name = "Gene set"

    class Input:
        """Input fields."""

        genes = ListField(
            StringField(),
            label="Genes",
            description="List of genes.",
        )
        source = StringField(
            label="Gene ID source",
            allow_custom_choice=True,
            choices=[
                ("AFFY", "AFFY"),
                ("DICTYBASE", "DICTYBASE"),
                ("ENSEMBL", "ENSEMBL"),
                ("NCBI", "NCBI"),
                ("UCSC", "UCSC"),
            ],
        )
        species = StringField(
            label="Species",
            description="Species latin name.",
            allow_custom_choice=True,
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
                ("Odocoileus virginianus texanus", "Odocoileus virginianus texanus"),
                ("Solanum tuberosum", "Solanum tuberosum"),
            ],
        )

    class Output:
        """Output fields."""

        geneset = FileField(label="Gene set")
        geneset_json = JsonField(label="Gene set (JSON)")
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")

    def run(self, inputs, outputs):
        """Run the analysis."""
        geneset = sorted(set(inputs.genes))
        if len(inputs.genes) != len(geneset):
            self.warning("Removed duplicated genes.")

        # Save geneset to file
        save_geneset_to_file(geneset, "geneset.tab.gz")
        outputs.geneset = "geneset.tab.gz"

        # Save geneset to json
        save_geneset_to_json(geneset, "geneset.json")
        outputs.geneset_json = "geneset.json"

        outputs.source = inputs.source
        outputs.species = inputs.species

        # Set the descriptor schema of this object to geneset
        ds = DescriptorSchema.get_latest(slug="geneset")
        self.data.descriptor_schema = ds.id


class CreateGenesetVenn(Process):
    """Create a gene set from a Venn diagram."""

    slug = "create-geneset-venn"
    name = "Gene set (create from Venn diagram)"
    process_type = "data:geneset:venn"
    version = "1.3.2"
    category = "Import"
    scheduling_class = SchedulingClass.INTERACTIVE
    persistence = Persistence.RAW
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {
                "image": "public.ecr.aws/genialis/resolwebio/base:ubuntu-22.04-14112023"
            }
        },
        "resources": {
            "cores": 1,
            "memory": 1024,
            "storage": 10,
        },
    }
    data_name = "Gene set (Venn)"

    class Input:
        """Input fields."""

        genes = ListField(
            StringField(),
            label="Genes",
            description="List of genes.",
        )
        source = StringField(
            label="Gene ID source",
            allow_custom_choice=True,
            choices=[
                ("AFFY", "AFFY"),
                ("DICTYBASE", "DICTYBASE"),
                ("ENSEMBL", "ENSEMBL"),
                ("NCBI", "NCBI"),
                ("UCSC", "UCSC"),
            ],
        )
        species = StringField(
            label="Species",
            description="Species latin name.",
            allow_custom_choice=True,
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
                ("Odocoileus virginianus texanus", "Odocoileus virginianus texanus"),
                ("Solanum tuberosum", "Solanum tuberosum"),
            ],
        )
        venn = FileField(
            label="Venn diagram",
            description="JSON file of Venn diagram.",
        )

    class Output:
        """Output fields."""

        geneset = FileField(label="Gene set")
        geneset_json = JsonField(label="Gene set (JSON)")
        source = StringField(label="Gene ID source")
        species = StringField(label="Species")
        venn = JsonField(label="Venn diagram")

    def run(self, inputs, outputs):
        """Run the analysis."""
        geneset = sorted(set(inputs.genes))
        if len(inputs.genes) != len(geneset):
            self.warning("Removed duplicated genes.")

        # Save geneset to file
        save_geneset_to_file(geneset, "geneset.tab.gz")
        outputs.geneset = "geneset.tab.gz"

        # Save geneset to json
        save_geneset_to_json(geneset, "geneset.json")
        outputs.geneset_json = "geneset.json"

        outputs.source = inputs.source
        outputs.species = inputs.species

        venn_file = inputs.venn.import_file(imported_format="extracted")
        outputs.venn = venn_file

        # Set the descriptor schema of this object to geneset
        ds = DescriptorSchema.get_latest(slug="geneset")
        self.data.descriptor_schema = ds.id
