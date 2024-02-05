"""Upload caches for Ensembl-VEP."""

import shutil
from pathlib import Path

from resolwe.process import DirField, FileField, Process, SchedulingClass, StringField


class ImportEnsemblVepCache(Process):
    """Import VEP cache directory."""

    slug = "upload-vep-cache"
    name = "Ensembl-VEP cache directory"
    process_type = "data:vep:cache"
    version = "1.1.0"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/dnaseq:6.3.1"},
        },
    }
    data_name = '{{ cache_file.file|default("?") }}'

    class Input:
        """Input fields to process Import VEP cache directory."""

        cache_file = FileField(label="Compressed cache directory")
        species = StringField(
            label="Species",
            description="Select a species name from the dropdown menu.",
            default="Homo sapiens",
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
            ],
        )

        build = StringField(
            label="Genome build",
        )
        release = StringField(label="Cache release")

    class Output:
        """Output fields to process Import VEP cache directory."""

        cache = DirField(label="Cache directory")
        species = StringField(label="Species")
        build = StringField(label="Build")
        release = StringField(label="Cache release")

    def run(self, inputs, outputs):
        """Run the analysis."""
        cache_dir = Path("Ensembl-vep_cache")
        if not cache_dir.exists():
            cache_dir.mkdir()

        cache_file = inputs.cache_file.import_file(imported_format="compressed")
        shutil.unpack_archive(cache_file, cache_dir)

        outputs.cache = cache_dir.name
        outputs.species = inputs.species
        outputs.build = inputs.build
        outputs.release = inputs.release
