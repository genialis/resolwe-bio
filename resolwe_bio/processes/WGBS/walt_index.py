"""Create genome index for WALT aligner."""
import shutil
from pathlib import Path

from plumbum import TEE

from resolwe.process import Cmd, DataField, DirField, FileField, Process, StringField


class WaltIndex(Process):
    """Create WALT genome index."""

    slug = "walt-index"
    process_type = "data:index:walt"
    name = "WALT genome index"
    requirements = {
        "expression-engine": "jinja",
        "executor": {"docker": {"image": "resolwebio/rnaseq:4.9.0"},},
        "resources": {"cores": 1, "memory": 16384,},
    }
    category = "Genome index"
    data_name = '{{ ref_seq.fasta.file|basename|default("?") }}'
    version = "1.0.1"

    class Input:
        """Input fields for WaltIndex."""

        ref_seq = DataField(
            "seq:nucleotide", label="Reference sequence (nucleotide FASTA)"
        )

    class Output:
        """Output fields to process WaltIndex."""

        index = DirField(label="WALT index")
        fastagz = FileField(label="FASTA file (compressed)")
        fasta = FileField(label="FASTA file")
        fai = FileField(label="FASTA file index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = Path(inputs.ref_seq.fasta.path).name
        assert basename.endswith(".fasta")
        name = basename[:-6]

        index_dir = Path("walt_index")
        index_dir.mkdir()

        shutil.copy(Path(inputs.ref_seq.fasta.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.fastagz.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.fai.path), Path.cwd())

        args = [
            "-c",
            inputs.ref_seq.fasta.path,
            "-o",
            index_dir / f"{name}.dbindex",
        ]

        return_code, _, _ = Cmd["makedb-walt"][args] & TEE(retcode=None)
        if return_code:
            self.error("Error occurred while preparing the WALT index.")

        outputs.index = index_dir.name
        outputs.fasta = f"{name}.fasta"
        outputs.fastagz = f"{name}.fasta.gz"
        outputs.fai = f"{name}.fasta.fai"
        outputs.species = inputs.ref_seq.species
        outputs.build = inputs.ref_seq.build
