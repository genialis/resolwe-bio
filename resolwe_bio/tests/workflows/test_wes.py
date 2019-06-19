# pylint: disable=missing-docstring
from resolwe.test import tag_process
from resolwe.flow.models import Data
from resolwe_bio.utils.test import BioProcessTestCase


class WESTestCase(BioProcessTestCase):
    @tag_process('workflow-wes')
    def test_wes(self):
        species = 'Homo sapiens'
        build = 'custom_build'

        with self.preparation_stage():
            reads = self.prepare_paired_reads(
                mate1=['./workflow_wes/input/STK11_1.fastq.gz'],
                mate2=['./workflow_wes/input/STK11_2.fastq.gz']
            )

            genome = self.run_process('upload-genome', {
                'src': './bqsr/input/hs_b37_chr19_upto_stk11.fasta.gz',
                'species': species,
                'build': build
            })

            bc_bedpe = self.run_process('upload-bedpe', {
                'src': './bamclipper/input/STK11.bedpe',
                'species': species,
                'build': build
            })

            kbase = []
            for i in ['./bqsr/input/dbsnp_STK11.vcf.gz']:
                kbase.append(
                    self.run_process('upload-variants-vcf', {
                        'src': i, 'species': species, 'build': build
                    }))

            intervals = self.run_process('upload-bed', {
                'src': './bqsr/input/STK11.bed',
                'species': species, 'build': build})

            adapters = self.prepare_adapters()

        input_workflow = {
            'reads': reads.pk,
            'genome': genome.id,
            'bamclipper_bedpe': bc_bedpe.id,
            'known_sites': [i.id for i in kbase],
            'intervals': intervals.id,
            'hc_dbsnp': kbase[0].id,
            'advanced': {
                'trimming': {
                    'adapters': adapters.id,
                    'seed_mismatches': 2,
                    'simple_clip_threshold': 10,
                    'min_adapter_length': 8,
                    'palindrome_clip_threshold': 30,
                    'leading': 20,
                    'trailing': 3,
                    'minlen': 40
                },
                'align': {
                    'seed_l': 19,
                    'band_w': 100,
                    'm': True,
                    're_seeding': 1.5,
                    'scoring': {
                        'match': 1,
                        'mismatch': 4,
                        'gap_o': 6,
                        'gap_e': 1,
                        'clipping': 5,
                        'unpaired_p': 17
                    },
                    'report_tr': 30
                },
                'markduplicates': {
                    'md_skip': False,
                    'md_remove_duplicates': False,
                    'md_validation_stringency': 'STRICT',
                    'md_assume_sort_order': ''
                },
                'bqsr': {
                    'read_group': '-LB=DAB;-PL=Illumina;-PU=barcode;-SM=sample1'
                },
                'hc': {
                    'stand_call_conf': 3,
                    'mbq': 3
                }
            }
        }

        self.run_process('workflow-wes', input_workflow)
        wes = Data.objects.last()
        self.assertFileExists(wes, 'vcf')
