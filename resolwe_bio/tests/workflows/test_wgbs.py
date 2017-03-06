# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase
from resolwe.flow.models import Data


class WgbsWorkflowTestCase(BioProcessTestCase):
    def test_wgbs_workflow(self):
        inputs = {'src': 'chr1_part.fasta.gz'}
        genome = self.run_process('upload-genome', inputs)

        inputs = {'src': ['wgbs.fastq.gz']}
        reads = self.run_process('upload-fastq-single', inputs)

        self.run_process(
            'workflow-wgbs', {
                'genome': genome.pk,
                'reads': reads.pk,
                'genome_identifier': 'hg19',
                'threads': 3
            }
        )

        wgbs = Data.objects.last()
        self.assertFile(wgbs, 'stats', 'wgbs.bam_stat.txt')
