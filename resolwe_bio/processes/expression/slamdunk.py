"""Run slamdunk all on paired-end reads."""
import os

from resolwe.process import Process, Cmd, DataField, BooleanField, IntegerField, FileField


class SlamdunkAllPaired(Process):
    """Run slamdunk for paired-end reads."""

    slug = 'slamdunk-all-paired'
    process_type = 'data:slamdunk:all'
    name = 'Slamdunk(paired-end)'
    requirements = {
        'expression-engine': 'jinja',
        'executor': {
            'docker': {
                'image': 'resolwebio/slamdunk:1.0.0'
            },
        },
        'resources': {
            'cores': 8,
            'memory': 16384,
        },
    }
    data_name = '{{ reads|sample_name|default("?") }}'
    version = '1.0.0'

    class Input:
        """Input fields for SlamdunkAllPaired."""

        reads = DataField('reads:fastq:paired', label='Reads')
        transcriptome = DataField('seq:nucleotide', label='FASTA file containig sequences for alingnig.')
        regions = DataField('bed', label='BED file with coordinates of regions of interest.')

        filter_multimappers = BooleanField(
            label='Filter multimappers',
            description='If true filter and reasign multimappers based on provided BED file with regions of interest.',
            default=True
        )

        max_alignments = IntegerField(
            label='Maximum number of multimapper alignments',
            description='The maximum number of alignments that will be reported for a multi-mapping read (i.e. reads'
                        'with multiple alignments of equal best scores).',
            default=1
        )

        read_length = IntegerField(
            label='Maximum read length',
            description='Maximul length of reads in the input FASTQ file.',
            default=150
        )

    class Output:
        """Output fields to process SlamdunkAllPaired."""

        tcount = FileField(label='Count report containing SLAMSeq statistics')
        map_output = FileField(label='Output of Slamdunk map')
        filter_output = FileField(label='Output of Slamdunk filter step')
        snp_output = FileField(label='Output of Slamdunk variant calls')
        count_output = FileField(label='Output of Slamdunk count')

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.reads.fastq[0].path)
        assert basename.endswith('fastq.gz')
        name = basename[:-9]

        reads_paths = [reads.path for reads in inputs.reads.fastq] + [reads.path for reads in inputs.reads.fastq2]
        (Cmd['ln']['-s', inputs.transcriptome.fasta.path, 'transcriptome.fasta'])()
        concatenated_name = f'{name}_concatenated.fastq.gz'
        (Cmd['cat'][reads_paths] > concatenated_name)()

        args = [
            '-r', 'transcriptome.fasta',
            '-b', inputs.regions.bed.path,
            '-o', 'slamdunk_all',
            '-n', inputs.max_alignments,
            '-t', '8',
            '-rl', inputs.read_length
        ]

        if inputs.filter_multimappers:
            args.extend(['-m', '-fb', inputs.regions.bed.path])

        args.append(concatenated_name)

        (Cmd['slamdunk']['all'][args])()

        Cmd['zip']['-r', f'{name}_map.zip', os.path.join('slamdunk_all', 'map')]()
        Cmd['zip']['-r', f'{name}_filter.zip', os.path.join('slamdunk_all', 'filter')]()
        Cmd['zip']['-r', f'{name}_snp.zip', os.path.join('slamdunk_all', 'snp')]()
        Cmd['zip']['-r', f'{name}_count.zip', os.path.join('slamdunk_all', 'count')]()

        outputs.tcount = os.path.join(
            'slamdunk_all',
            'count',
            f'{name}_concatenated.fastq_slamdunk_mapped_filtered_tcount.tsv'
        )
        outputs.map_output = f'{name}_map.zip'
        outputs.filter_output = f'{name}_filter.zip'
        outputs.snp_output = f'{name}_snp.zip'
        outputs.count_output = f'{name}_count.zip'
