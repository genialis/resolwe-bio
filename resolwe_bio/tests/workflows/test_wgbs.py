# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase
from resolwe_bio.utils.filter import filter_comment_lines


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

        hmr = Data.objects.filter(process__slug='hmr').last()
        self.assertFile(hmr, 'hmr', '3A_WT_WGBS_chr2_17k_single.hmr.gz', compression='gzip')
        summary = Data.objects.filter(process__slug='alignment-summary').last()
        self.assertFile(
            summary,
            'report',
            '3A_WT_WGBS_alignment_metrics.txt',
            file_filter=filter_comment_lines)
        wgs_metrics = Data.objects.filter(process__slug='wgs-metrics').last()
        self.assertFile(
            wgs_metrics,
            'report',
            '3A_WT_WGBS_wgs_metrics.txt',
            file_filter=filter_comment_lines)
        rrbs_metrics = Data.objects.filter(process__slug='rrbs-metrics').last()
        self.assertFile(
            rrbs_metrics,
            'report',
            '3A_WT_WGBS_rrbs_summary_metrics.txt',
            file_filter=filter_comment_lines)
