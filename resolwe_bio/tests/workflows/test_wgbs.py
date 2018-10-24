# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class WgbsWorkflowTestCase(BioProcessTestCase):
    @tag_process('workflow-wgbs')
    def test_wgbs_workflow(self):
        with self.preparation_stage():
            inputs = {
                'src': 'hg19_chr2 17k.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            genome = self.run_process('upload-genome', inputs)

            inputs = {'src': ['3A_WT_WGBS_chr2_17k R1.fastq.gz']}
            reads = self.run_process('upload-fastq-single', inputs)

        self.run_process(
            'workflow-wgbs', {
                'genome': genome.id,
                'reads': reads.id,
                'alignment': {'rm_dup': False},
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        wgbs = Data.objects.last()
        self.assertFile(wgbs, 'hmr', '3A_WT_WGBS_chr2_17k_single.hmr.gz', compression='gzip')
