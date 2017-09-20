# pylint: disable=missing-docstring
from resolwe_bio.utils.filter import filter_vcf_variable
from resolwe_bio.utils.test import BioProcessTestCase
from resolwe.flow.models import Data


class CheMutWorkflowTestCase(BioProcessTestCase):

    def test_chemut_workflow(self):
        inputs = {'src': 'chemut_genome.fasta.gz'}
        genome = self.run_process('upload-genome', inputs)

        inputs = {'src1': ['AX4_mate1.fq.gz'],
                  'src2': ['AX4_mate2.fq.gz']}

        parental_reads = self.run_process('upload-fastq-paired', inputs)

        inputs = {'src1': ['CM_mate1.fq.gz'],
                  'src2': ['CM_mate2.fq.gz']}

        mut_reads = self.run_process('upload-fastq-paired', inputs)

        inputs = {'genome': genome.id, 'reads': parental_reads.id}
        align_parental = self.run_process('alignment-bwa-mem', inputs)

        inputs = {'genome': genome.id, 'reads': mut_reads.id}
        align_mut = self.run_process('alignment-bwa-mem', inputs)

        self.run_process(
            'workflow-chemut', {
                'analysis_type': 'snv',
                'parental_strains': [align_parental.id],
                'mutant_strains': [align_mut.id],
                'genome': genome.id
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        variants = Data.objects.last()
        self.assertFile(variants, 'vcf', 'chemut.vcf', file_filter=filter_vcf_variable)
