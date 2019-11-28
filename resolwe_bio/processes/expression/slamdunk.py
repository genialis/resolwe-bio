"""Run slamdunk all on paired-end reads."""
import glob
import os

from plumbum import TEE

from resolwe.process import BooleanField, Cmd, DataField, DirField, IntegerField, FileField, Process, StringField


class SlamdunkAllPaired(Process):
    """Run slamdunk for paired-end reads."""

    slug = 'slamdunk-all-paired'
    process_type = 'data:alignment:bam:slamdunk'
    name = 'Slamdunk (paired-end)'
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
    entity = {
        'type': 'sample',
    }
    data_name = '{{ reads|sample_name|default("?") }}'
    version = '2.0.0'

    class Input:
        """Input fields for SlamdunkAllPaired."""

        reads = DataField('reads:fastq:paired', label='Reads')
        ref_seq = DataField('seq:nucleotide', label='FASTA file.')
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
        bam = FileField(label='Aligned reads (BAM)')
        bai = FileField(label='Aligned reads index (BAI)')
        bigwig = FileField(label="BigWig file", required=False)
        stats = FileField(label='Alignment statistics')
        map_output = FileField(label='Output of Slamdunk map step')
        map_output_dir = DirField(label='Output of Slamdunk map step')
        filter_output = FileField(label='Output of Slamdunk filter step')
        filter_output_dir = DirField(label='Output of Slamdunk filter step')
        snp_output = FileField(label='Output of Slamdunk variant calls')
        snp_output_dir = DirField(label='Output of Slamdunk variant calls')
        count_output = FileField(label='Output of Slamdunk count')
        count_output_dir = DirField(label='Output of Slamdunk count')
        species = StringField(label="Species")
        build = StringField(label="Build")

    def run(self, inputs, outputs):
        """Run analysis."""
        basename = os.path.basename(inputs.reads.fastq[0].path)
        assert basename.endswith('fastq.gz')
        name = basename[:-9]

        reads_paths = [reads.path for reads in inputs.reads.fastq] + [reads.path for reads in inputs.reads.fastq2]
        (Cmd['ln']['-s', inputs.ref_seq.fasta.path, 'ref_seq.fasta'])()
        concatenated_name = f'{name}_concatenated.fastq.gz'
        (Cmd['cat'][reads_paths] > concatenated_name)()

        args = [
            '-r', 'ref_seq.fasta',
            '-b', inputs.regions.bed.path,
            '-o', 'slamdunk_all',
            '-n', inputs.max_alignments,
            '-t', self.requirements.resources.cores,
            '-rl', inputs.read_length,
        ]

        if inputs.filter_multimappers:
            args.extend(['-m', '-fb', inputs.regions.bed.path])

        args.append(concatenated_name)

        return_code, _, _ = Cmd['slamdunk']['all'][args] & TEE(retcode=None)
        if return_code:
            self.error('Slamdunk analysis failed.')

        bam_path = os.path.join('slamdunk_all', 'map', '*.bam')
        bam_file = glob.glob(bam_path)[0]
        if os.path.isfile(bam_file):
            bam_name = os.path.basename(bam_file)[:-4]
            sorted_bam = bam_name + '_sorted.bam'
            (Cmd['samtools']['sort', bam_file] > sorted_bam)()
            Cmd['samtools']['index', sorted_bam]()
            (Cmd['samtools']['flagstat', sorted_bam] > bam_name + '_sorted_stats.txt')()
            Cmd['bamtobigwig.sh']([sorted_bam, inputs.regions.species, self.requirements.resources.cores])
            outputs.bam = sorted_bam
            outputs.bai = sorted_bam + '.bai'
            outputs.stats = bam_name + '_sorted_stats.txt'
            outputs.bigwig = bam_name + '_sorted.bw'
        else:
            self.error('Slamdunk failed to create a BAM file. Exiting.')

        tcount_path = os.path.join('slamdunk_all', 'count', '*_tcount.tsv')
        tcount_file = glob.glob(tcount_path)[0]
        if os.path.isfile(tcount_file):
            outputs.tcount = tcount_file
        else:
            self.error('Slamdunk failed to create a tcount file. Exiting.')

        outputs.map_output_dir = os.path.join('slamdunk_all', 'map')
        outputs.filter_output_dir = os.path.join('slamdunk_all', 'filter')
        outputs.snp_output_dir = os.path.join('slamdunk_all', 'snp')
        outputs.count_output_dir = os.path.join('slamdunk_all', 'count')

        Cmd['zip']['-r', f'{name}_map.zip', os.path.join('slamdunk_all', 'map')]()
        Cmd['zip']['-r', f'{name}_filter.zip', os.path.join('slamdunk_all', 'filter')]()
        Cmd['zip']['-r', f'{name}_snp.zip', os.path.join('slamdunk_all', 'snp')]()
        Cmd['zip']['-r', f'{name}_count.zip', os.path.join('slamdunk_all', 'count')]()

        outputs.map_output = f'{name}_map.zip'
        outputs.filter_output = f'{name}_filter.zip'
        outputs.snp_output = f'{name}_snp.zip'
        outputs.count_output = f'{name}_count.zip'
        outputs.species = inputs.regions.species
        outputs.build = inputs.regions.build
