"""Import Nucleotide Sequence (FASTA)."""
from pathlib import Path

import dnaio
from dnaio.exceptions import FastaFormatError

from resolwe.process import (
    Cmd,
    FileField,
    IntegerField,
    Persistence,
    Process,
    SchedulingClass,
    StringField,
)


class ImportFastaNucleotide(Process):
    """Import nucleotide sequence file in FASTA format.

    FASTA file is a text-based format for representing nucleotide sequences, in which nucleotides
    are represented using single-letter codes. The uploaded FASTA file can hold multiple nucleotide
    sequences.
    """

    slug = "upload-fasta-nucl"
    name = "FASTA file"
    process_type = "data:seq:nucleotide"
    version = "3.0.1"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    persistence = Persistence.RAW
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/rnaseq:4.9.0"}},
        "resources": {
            "cores": 2,
            "memory": 8192,
            "network": True,
        },
    }
    data_name = '{{ src.file|default("?") }}'

    class Input:
        """Input fields to process ImportFastaNucleotide."""

        src = FileField(label="Sequence file (FASTA)")

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
            description="Enter a genome build information associated "
            "with the uploaded sequence(s).",
        )

    class Output:
        """Output field of the process ImportFastaNucleotide."""

        fastagz = FileField(label="FASTA file (compressed)")
        fasta = FileField(label="FASTA file")
        fai = FileField(label="FASTA file index")
        fasta_dict = FileField(label="FASTA dictionary")
        num_seqs = IntegerField(label="Number of sequences")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def validate_fasta(self, infile):
        """Validate FASTA file format."""
        try:
            with dnaio.open(infile, fileformat="fasta") as fasta_file:
                if not any(fasta_file):
                    self.error(
                        f"The uploaded .FASTA file {infile} contains no sequence data."
                    )
                else:
                    self.info(f"Successfully validated the input file {infile}.")
        except FastaFormatError as dnaio_error:
            self.error(f"Format error in the uploaded file {infile}. {dnaio_error}")

    def run(self, inputs, outputs):
        """Run the analysis."""
        fasta_path = inputs.src.import_file(
            imported_format="extracted", progress_to=0.2
        )

        supported_extensions = (".fa", ".fasta", ".faa", ".fna", ".ffn", ".frn")
        if not fasta_path.endswith(supported_extensions):
            self.error(
                f"The imported file has unsupported file name extension. "
                f"The supported extensions are {supported_extensions}."
            )

        # validate the format of the uploaded fasta file
        self.validate_fasta(fasta_path)
        self.progress = 0.2

        # ensure the .fasta suffix and compress the output
        fasta_path = Path(fasta_path)
        output_fasta = fasta_path.with_suffix(".fasta")
        fasta_path.rename(output_fasta)
        Cmd["pigz"]["-k", output_fasta]()
        self.progress = 0.4

        # create a .fai index file
        Cmd["samtools"]["faidx", output_fasta]()
        self.progress = 0.6

        # Create fasta dictionary file
        fasta_dict = f"{output_fasta.stem}.dict"
        Cmd["picard-tools"][
            "CreateSequenceDictionary", f"R={output_fasta.name}", f"O={fasta_dict}"
        ]()
        self.progress = 0.8

        # check the number of sequences in the .fasta file
        seq_number = Cmd["grep"]["-c", "^>", output_fasta]().strip()

        # save the outputs
        outputs.fasta = output_fasta.name
        outputs.fastagz = f"{output_fasta}.gz"
        outputs.fai = f"{output_fasta}.fai"
        outputs.fasta_dict = fasta_dict
        outputs.num_seqs = seq_number
        outputs.species = inputs.species
        outputs.build = inputs.build
