"""Create genome index for BWA aligner."""

import shutil
from pathlib import Path

from plumbum import TEE

from resolwe.process import (
    Cmd,
    DataField,
    DirField,
    FileField,
    Persistence,
    Process,
    StringField,
)


class BWAIndex(Process):
    """Create BWA genome index."""

    slug = "bwa-index"
    process_type = "data:index:bwa"
    name = "BWA genome index"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/rnaseq:6.0.0"},
        },
        "resources": {
            "cores": 1,
            "memory": 16384,
        },
    }
    category = "Genome index"
    data_name = '{{ ref_seq.fasta.file|basename|default("?") }}'
    version = "1.2.1"
    persistence = Persistence.CACHED

    class Input:
        """Input fields for BWAIndex."""

        ref_seq = DataField(
            "seq:nucleotide", label="Reference sequence (nucleotide FASTA)"
        )

    class Output:
        """Output fields to process BWAIndex."""

        index = DirField(label="BWA index")
        fastagz = FileField(label="FASTA file (compressed)")
        fasta = FileField(label="FASTA file")
        fai = FileField(label="FASTA file index")
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = Path(inputs.ref_seq.output.fasta.path).name
        assert basename.endswith(".fasta")
        name = basename[:-6]

        index_dir = Path("BWA_index")
        index_dir.mkdir()

        shutil.copy(Path(inputs.ref_seq.output.fasta.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.output.fastagz.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.output.fai.path), Path.cwd())

        args = [
            "-p",
            index_dir / f"{name}.fasta",
            inputs.ref_seq.output.fasta.path,
        ]

        return_code, _, _ = Cmd["bwa"]["index"][args] & TEE(retcode=None)
        if return_code:
            self.error("Error occurred while preparing the BWA index.")

        outputs.index = index_dir.name
        outputs.fasta = f"{name}.fasta"
        outputs.fastagz = f"{name}.fasta.gz"
        outputs.fai = f"{name}.fasta.fai"
        outputs.species = inputs.ref_seq.output.species
        outputs.build = inputs.ref_seq.output.build
