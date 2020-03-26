"""Create genome index for Subread aligner."""
import shutil

from pathlib import Path
from plumbum import TEE
from resolwe.process import Cmd, DataField, FileField, Process, StringField, DirField


class SubreadIndex(Process):
    """Create Subread genome index."""

    slug = 'subread-index'
    process_type = 'data:index:subread'
    name = 'Subread genome index'
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/rnaseq:4.9.0'
            },
        },
        'resources': {
            'cores': 1,
            'memory': 16384,
        },
    }
    category = 'Genome index'
    data_name = '{{ ref_seq.file|default("?") }}'
    version = '1.0.0'

    class Input:
        """Input fields for SubreadIndex."""

        ref_seq = DataField('seq:nucleotide', label='Reference sequence (nucleotide FASTA)')

    class Output:
        """Output fields to process SubreadIndex."""

        index = DirField(label='Subread index')
        fastagz = FileField(label='FASTA file (compressed)')
        fasta = FileField(label='FASTA file')
        fai = FileField(label='FASTA file index')
        species = StringField(label='Species')
        build = StringField(label='Build')

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = Path(inputs.ref_seq.fasta.path).name
        assert basename.endswith('.fasta')
        name = basename[:-6]

        index_dir = Path('subread_index')
        index_dir.mkdir()

        shutil.copy(Path(inputs.ref_seq.fasta.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.fastagz.path), Path.cwd())
        shutil.copy(Path(inputs.ref_seq.fai.path), Path.cwd())

        args = [
            inputs.ref_seq.fasta.path,
            '-o', index_dir / f'{name}_index',
        ]

        return_code, _, _ = Cmd['subread-buildindex'][args] & TEE(retcode=None)
        if return_code:
            self.error("Error occurred while preparing the Subread index.")

        outputs.index = index_dir.name
        outputs.fasta = f'{name}.fasta'
        outputs.fastagz = f'{name}.fasta.gz'
        outputs.fai = f'{name}.fasta.fai'
        outputs.species = inputs.ref_seq.species
        outputs.build = inputs.ref_seq.build
