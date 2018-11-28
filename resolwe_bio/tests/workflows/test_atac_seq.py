# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class AtacSeqWorkflowTestCase(BioProcessTestCase):
    @tag_process('workflow-atac-seq')
    def test_atac_seq_workflow(self):
        with self.preparation_stage():
            inputs = {
                'src': 'mm10_chr17_44M-45M.fa.gz',
                'species': 'Mus musculus',
                'build': 'mm10'
            }
            genome = self.run_process('upload-genome', inputs)
            reads = self.prepare_paired_reads(mate1=['atac_R1.fastq.gz'],
                                              mate2=['atac_R2.fastq.gz'])
            inputs = {
                'src': 'promoters_mm10_chr17_subregion.bed',
                'species': 'Mus musculus',
                'build': 'mm10',
            }
            promoters = self.run_process('upload-bed', inputs)

        self.run_process(
            'workflow-atac-seq', {
                'reads': reads.id,
                'genome': genome.id,
                'promoter': promoters.id,
                'settings': {'bedgraph': False},
            }
        )
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        macs2 = Data.objects.last()
        self.assertFile(macs2, 'chip_qc', 'atac_seq_report.txt')
