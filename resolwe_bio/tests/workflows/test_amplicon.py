# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class AmpliconWorkflowTestCase(BioProcessTestCase):
    @tag_process('workflow-accel')
    def test_amplicon_workflow(self):
        with self.preparation_stage():
            inputs = {
                'src1': ['56GSID_10k_mate1.fastq.gz'],
                'src2': ['56GSID_10k_mate2.fastq.gz']}
            reads = self.run_process('upload-fastq-paired', inputs)

            adapters = self.run_process('upload-fasta-nucl', {'src': 'adapters.fasta'})

            inputs = {
                'src': 'hs_b37_chr2_small.fasta.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            genome = self.run_process('upload-genome', inputs)

            master_file = self.prepare_amplicon_master_file()

            vcf_input = {
                'src': '1000G_phase1.indels.b37_chr2_small.vcf.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            indels = self.run_process('upload-variants-vcf', vcf_input)

            dbsnp_input = {
                'src': 'dbsnp_138.b37.chr2_small.vcf.gz',
                'species': 'Homo sapiens',
                'build': 'b37'
            }
            dbsnp = self.run_process('upload-variants-vcf', dbsnp_input)

        self.run_process(
            'workflow-accel', {
                'reads': reads.id,
                'genome': genome.id,
                'master_file': master_file.id,
                'adapters': adapters.id,
                'preprocess_bam': {
                    'known_vars': [dbsnp.id],
                    'known_indels': [indels.id]
                },
                'gatk': {
                    'dbsnp': dbsnp.id,
                    'mbq': 20,
                    'stand_call_conf': 20
                },
                'lofreq': {
                    'min_bq': 20,
                    'min_alt_bq': 20
                },
                'var_annot': {
                    'known_vars_db': [dbsnp.id]
                },
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
