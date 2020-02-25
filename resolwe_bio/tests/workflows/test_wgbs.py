# pylint: disable=missing-docstring
import os

from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase
from resolwe_bio.utils.filter import filter_comment_lines


class WgbsWorkflowTestCase(BioProcessTestCase):
    @tag_process('workflow-wgbs-single', 'workflow-wgbs-paired')
    def test_wgbs_workflow(self):
        with self.preparation_stage():
            inputs = {
                'src': os.path.join('wgbs', 'input', 'hg19_chr2 17k.fasta.gz'),
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            genome = self.run_process('upload-genome', inputs)

            inputs = {'src': [os.path.join('wgbs', 'input', '3A_WT_WGBS_chr2_17k_R1.fastq.gz')]}
            reads = self.run_process('upload-fastq-single', inputs)
            reads_paired = self.prepare_paired_reads(
                mate1=[os.path.join('wgbs', 'input', '3A_WT_WGBS_chr2_17k_R1.fastq.gz')],
                mate2=[os.path.join('wgbs', 'input', '3A_WT_WGBS_chr2_17k_R2.fastq.gz')]
            )

        self.run_process(
            'workflow-wgbs-single', {
                'genome': genome.id,
                'reads': reads.id,
                'alignment': {'rm_dup': False},
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        hmr = Data.objects.filter(process__slug='hmr').last()
        self.assertFile(
            hmr,
            'hmr',
            os.path.join('wgbs_workflow', 'output', 'single_end', '3A_WT_WGBS_chr2_17k_single.hmr.gz'),
            compression='gzip'
        )
        summary = Data.objects.filter(process__slug='alignment-summary').last()
        self.assertFile(
            summary,
            'report',
            os.path.join('wgbs_workflow', 'output', 'single_end', '3A_WT_WGBS_alignment_metrics.txt'),
            file_filter=filter_comment_lines
        )
        wgs_metrics = Data.objects.filter(process__slug='wgs-metrics').last()
        self.assertFile(
            wgs_metrics,
            'report',
            os.path.join('wgbs_workflow', 'output', 'single_end', '3A_WT_WGBS_wgs_metrics.txt'),
            file_filter=filter_comment_lines
        )
        rrbs_metrics = Data.objects.filter(process__slug='rrbs-metrics').last()
        self.assertFile(
            rrbs_metrics,
            'report',
            os.path.join('wgbs_workflow', 'output', 'single_end', '3A_WT_WGBS_rrbs_summary_metrics.txt'),
            file_filter=filter_comment_lines
        )

        self.run_process(
            'workflow-wgbs-paired', {
                'genome': genome.id,
                'reads': reads_paired.id,
                'alignment': {'rm_dup': False},
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        hmr = Data.objects.filter(process__slug='hmr').last()
        self.assertFile(
            hmr,
            'hmr',
            os.path.join('wgbs_workflow', 'output', 'paired_end', '3A_WT_WGBS_chr2_17k_paired.hmr.gz'),
            compression='gzip'
        )
        summary = Data.objects.filter(process__slug='alignment-summary').last()
        self.assertFile(
            summary,
            'report',
            os.path.join('wgbs_workflow', 'output', 'paired_end', '3A_WT_WGBS_paired_alignment_metrics.txt'),
            file_filter=filter_comment_lines
        )
        wgs_metrics = Data.objects.filter(process__slug='wgs-metrics').last()
        self.assertFile(
            wgs_metrics,
            'report',
            os.path.join('wgbs_workflow', 'output', 'paired_end', '3A_WT_WGBS_paired_wgs_metrics.txt'),
            file_filter=filter_comment_lines
        )
        rrbs_metrics = Data.objects.filter(process__slug='rrbs-metrics').last()
        self.assertFile(
            rrbs_metrics,
            'report',
            os.path.join('wgbs_workflow', 'output', 'paired_end', '3A_WT_WGBS_paired_rrbs_summary_metrics.txt'),
            file_filter=filter_comment_lines
        )
        insert_size = Data.objects.filter(process__slug='insert-size').last()
        self.assertFile(
            insert_size,
            'report',
            os.path.join('wgbs_workflow', 'output', 'paired_end', '3A_WT_WGBS_paired_insert_size_metrics.txt'),
            file_filter=filter_comment_lines
        )
