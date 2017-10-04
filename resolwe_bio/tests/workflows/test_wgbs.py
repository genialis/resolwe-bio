# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class WgbsWorkflowTestCase(BioProcessTestCase):
    @tag_process('workflow-wgbs')
    def test_wgbs_workflow(self):
        with self.preparation_stage():
            inputs = {
                'src': 'chr1_part.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            genome = self.run_process('upload-genome', inputs)

            inputs = {'src': ['wgbs.fastq.gz']}
            reads = self.run_process('upload-fastq-single', inputs)

        self.run_process(
            'workflow-wgbs', {
                'genome': genome.id,
                'reads': reads.id,
                'genome_identifier': 'hg19',
                'threads': 3
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        wgbs = Data.objects.last()
        self.assertFile(wgbs, 'stats', 'wgbs.bam_stat.txt')
