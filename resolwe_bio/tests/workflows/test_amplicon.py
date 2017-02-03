# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class AmpliconWorkflowTestCase(BioProcessTestCase):
    @skipDockerFailure("Processor requires a custom Docker image.")
    def test_amplicon_workflow(self):
        inputs = {
            'src1': ['56GSID_10k_mate1.fastq.gz'],
            'src2': ['56GSID_10k_mate2.fastq.gz']}
        reads = self.run_process('upload-fastq-paired', inputs)

        adapters = self.run_process('upload-fasta-nucl', {'src': 'adapters.fasta'})
        genome = self.run_process('upload-genome', {'src': 'hs_b37_chr22_frag.fasta.gz'})

        inputs = {'src': '5ptrim_new56Gprimers.fa.gz'}
        primers_1 = self.run_processor('upload-fasta-nucl', inputs)

        inputs = {'src': '3ptrim_new56Gprimers.fa.gz'}
        primers_2 = self.run_processor('upload-fasta-nucl', inputs)

        intervals_picard = self.run_process('upload-bed', {'src': '56g_targets_picard_small.bed'})
        intervals = self.run_process('upload-bed', {'src': '56g_targets_small.bed'})

        inputs = {'src': 'Mills_and_1000G_gold_standard.indels.b37.chr22_small.vcf.gz'}
        indels_1 = self.run_process('upload-variants-vcf', inputs)

        inputs = {'src': '1000G_phase1.indels.b37.chr22_small.vcf.gz'}
        indels_2 = self.run_process('upload-variants-vcf', inputs)

        dbsnp = self.run_process('upload-variants-vcf', {'src': 'dbsnp_138.b37.chr22_small.vcf.gz'})

        self.run_process(
            'workflow-accel', {
                'reads': reads.id,
                'genome': genome.id,
                'primers': {
                    'adapters': adapters.id,
                    'up_primers': primers_1.id,
                    'down_primers': primers_2.id,
                },
                'target_intervals': {
                    'intervals_picard': intervals_picard.id,
                    'intervals': intervals.id,
                },
                'preprocess_bam': {
                    'known_vars': [dbsnp.id],
                    'known_indels': [indels_1.id, indels_2.id]
                },
                'gatk': {
                    'dbsnp': dbsnp.id
                },
                'var_annot': {
                    'known_vars_db': [dbsnp.id]
                }
            }
        )
