# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.filter import filter_vcf_variable
from resolwe_bio.utils.test import BioProcessTestCase


class HeatSeqWorkflowTestCase(BioProcessTestCase):
    @tag_process('workflow-heat-seq')
    def test_heatseq_workflow(self):
        with self.preparation_stage():
            inputs = {
                'src1': ['heat_seq_mate1.fq.gz'],
                'src2': ['heat_seq_mate2.fq.gz']}
            reads = self.run_process('upload-fastq-paired', inputs)

            inputs = {
                'src': 'chr21_small.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            genome = self.run_process('upload-genome', inputs)

            probe = self.run_process('upload-file', {'src': 'heat_seq_probe_info.txt'})
            bed_input = {
                'src': 'heat_seq_capture_targets.bed',
                'species': 'Homo sapiens',
                'build': 'hg19'
            }
            bed = self.run_process('upload-bed', bed_input)

        self.run_process(
            'workflow-heat-seq', {
                'reads': reads.id,
                'genome': genome.id,
                'probe_info': probe.id,
                'bed': bed.id
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        variants = Data.objects.last()

        self.assertFile(
            variants,
            'vcf',
            'heat-seq.vcf.gz',
            file_filter=filter_vcf_variable,
            compression='gzip'
        )
