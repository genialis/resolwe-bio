"""Create genome index for BWA-MEM2 aligner."""

import shutil
from pathlib import Path

from plumbum import TEE

from resolwe.process import Cmd, DataField, DirField, FileField, Process, StringField


class BWAMEM2Index(Process):
    """Create BWA-MEM2 genome index."""

    slug = "bwamem2-index"
    process_type = "data:index:bwamem2"
    name = "BWA-MEM2 genome index"
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/genialis/resolwebio/dnaseq:6.3.1"},
        },
        "resources": {
            "cores": 6,
            "memory": 98304,
        },
    }
    category = "Genome index"
    data_name = '{{ ref_seq.fasta.file|basename|default("?") }}'
    version = "1.1.0"

    class Input:
        """Input fields for BWAMEM2Index."""

        ref_seq = DataField(
            "seq:nucleotide", label="Reference sequence (nucleotide FASTA)"
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
        """Run analysis."""
        basename = Path(inputs.ref_seq.output.fasta.path).name
        assert basename.endswith(".fasta")
        name = basename[:-6]

        index_dir = Path("BWAMEM2_index")
        index_dir.mkdir()

        shutil.copy(Path(inputs.ref_seq.output.fasta.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.output.fastagz.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.output.fai.path), Path.cwd())

        args = [
            "-p",
            index_dir / f"{name}.fasta",
            inputs.ref_seq.output.fasta.path,
        ]

        return_code, _, _ = Cmd["bwa-mem2"]["index"][args] & TEE(retcode=None)
        if return_code:
            self.error("Error occurred while preparing the BWA-MEM2 index.")

        outputs.index = index_dir.name
        outputs.fasta = f"{name}.fasta"
        outputs.fastagz = f"{name}.fasta.gz"
        outputs.fai = f"{name}.fasta.fai"
        outputs.species = inputs.ref_seq.output.species
        outputs.build = inputs.ref_seq.output.build
