"""Create genome index for HISAT2 aligner."""
import shutil
from pathlib import Path

from plumbum import TEE

from resolwe.process import Cmd, DataField, DirField, FileField, Process, StringField


class Hisat2Index(Process):
    """Create HISAT2 genome index."""

    slug = "hisat2-index"
    process_type = "data:index:hisat2"
    name = "HISAT2 genome index"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/rnaseq:5.9.0"},
        },
        "resources": {
            "cores": 10,
            "memory": 16384,
        },
    }
    category = "Genome index"
    data_name = '{{ ref_seq.fasta.file|basename|default("?") }}'
    version = "1.1.1"

    class Input:
        """Input fields for Hisat2Index."""

        ref_seq = DataField(
            "seq:nucleotide", label="Reference sequence (nucleotide FASTA)"
        )

    class Output:
        """Output fields to process Hisat2Index."""

        index = DirField(label="HISAT2 index")
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

        index_dir = Path("hisat2_index")
        index_dir.mkdir()

        shutil.copy(Path(inputs.ref_seq.output.fasta.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.output.fastagz.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.output.fai.path), Path.cwd())

        args = [
            inputs.ref_seq.output.fasta.path,
            index_dir / f"{name}_index",
            "-p",
            self.requirements.resources.cores,
        ]

        return_code, _, _ = Cmd["hisat2-build"][args] & TEE(retcode=None)
        if return_code:
            self.error("Error occurred while preparing the HISAT2 index.")

        outputs.index = index_dir.name
        outputs.fasta = f"{name}.fasta"
        outputs.fastagz = f"{name}.fasta.gz"
        outputs.fai = f"{name}.fasta.fai"
        outputs.species = inputs.ref_seq.output.species
        outputs.build = inputs.ref_seq.output.build
