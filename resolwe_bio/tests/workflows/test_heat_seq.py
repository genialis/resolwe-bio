# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase
from resolwe.flow.models import Data


class HeatSeqWorkflowTestCase(BioProcessTestCase):
    def test_heatseq_workflow(self):
        inputs = {
            'src1': ['heat_seq_mate1.fq.gz'],
            'src2': ['heat_seq_mate2.fq.gz']}
        reads = self.run_process('upload-fastq-paired', inputs)

        genome = self.run_process('upload-genome', {'src': 'chr21_small.fasta.gz'})
        probe = self.run_process('upload-file', {'src': 'heat_seq_probe_info.txt'})
        bed = self.run_process('upload-bed', {'src': 'heat_seq_capture_targets.bed'})

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

        def filter_version(line):
            if line.startswith(b"##samtoolsVersion") or line.startswith(b"##reference") or b"/data_all/" in line:
                return True

        self.assertFile(variants, 'vcf', 'heat-seq.vcf', file_filter=filter_version)
