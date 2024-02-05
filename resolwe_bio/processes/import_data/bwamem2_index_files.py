"""Upload index files for BWA-MEM2 aligner."""

import shutil
from pathlib import Path

from resolwe.process import (
    Cmd,
    DirField,
    FileField,
    Process,
    SchedulingClass,
    StringField,
)


class ImportBWA2Index(Process):
    """Import BWA-MEM2 index files."""

    slug = "upload-bwamem2-index"
    name = "BWA-MEM2 index files"
    process_type = "data:index:bwamem2"
    version = "1.1.0"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/dnaseq:6.3.1"},
        },
    }
    data_name = '{{ ref_seq.file|default("?") }}'
    version = "1.0.0"

    class Input:
        """Input fields to process Import BWA-MEM2 index."""

        ref_seq = FileField(label="Reference sequence (nucleotide FASTA)")
        index_name = FileField(label="BWA-MEM2 index files")
        species = StringField(
            label="Species",
            description="Select a species name from the dropdown menu "
            "or write a custom species name in the species "
            "field. For sequences that are not related to "
            "any particular species (e.g. adapters file), "
            "you can select the value Other.",
            allow_custom_choice=True,
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
                ("Rattus norvegicus", "Rattus norvegicus"),
                ("Macaca mulatta", "Macaca mulatta"),
                ("Dictyostelium discoideum", "Dictyostelium discoideum"),
                ("Other", "Other"),
            ],
        )

        build = StringField(
            label="Genome build",
        )

    class Output:
        """Output fields to process BWAMEM2Index."""

        index = DirField(label="BWA-MEM2 index")
        fastagz = FileField(label="FASTA file (compressed)")
        fasta = FileField(label="FASTA file")
        fai = FileField(label="FASTA file index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run the analysis."""
        index_dir = Path("BWAMEM2_index")
        if not index_dir.exists():
            index_dir.mkdir()

        index_file = inputs.index_name.import_file(imported_format="compressed")
        shutil.unpack_archive(index_file, index_dir)

        fasta_path = inputs.ref_seq.import_file(imported_format="extracted")

        supported_extensions = (".fa", ".fasta", ".faa", ".fna", ".ffn", ".frn")
        if not fasta_path.endswith(supported_extensions):
            self.error(
                f"The imported file has unsupported file name extension. "
                f"The supported extensions are {supported_extensions}."
            )

        fasta_path = Path(fasta_path)
        output_fasta = fasta_path.with_suffix(".fasta")
        fasta_path.rename(output_fasta)

        Cmd["pigz"]["-k", output_fasta]()
        Cmd["samtools"]["faidx", output_fasta]()

        outputs.index = index_dir.name
        outputs.fasta = output_fasta.name
        outputs.fastagz = f"{output_fasta}.gz"
        outputs.fai = f"{output_fasta}.fai"
        outputs.species = inputs.species
        outputs.build = inputs.build
