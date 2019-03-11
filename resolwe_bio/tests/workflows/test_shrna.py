# pylint: disable=missing-docstring,invalid-name
from resolwe.test import tag_process
from resolwe.flow.models import Data
from resolwe_bio.utils.test import KBBioProcessTestCase


class SHRNATestCase(KBBioProcessTestCase):
    @tag_process('workflow-trim-align-quant')
    def test_run_trim_align_quant(self):
        with self.preparation_stage():
            pf_in = './shrna_diffexp/input/'
            pf_out = './shrna_diffexp/output/'

            species = 'Homo sapiens'
            build = 'customshRNA'
            reads = self.prepare_reads([pf_in + 'SM31_ss.fastq.gz'])
            # Larger data ran on this reduced genome yield similar alignment statistics (~2 % aligned reads).
            genome = self.run_process('upload-genome', {'src': pf_in + 'SM31_library.fasta.gz',
                                                        'species': species,
                                                        'build': build})

        input_workflow = {
            'reads': reads.id,
            'trimming_options': {
                'up_primers_seq': ['TAGTGAAGCCACAGAT'],
                'down_primers_seq': ['TACTGCCTCGGA'],
                'error_rate_5end': 0.2,
                'error_rate_3end': 0.2
            },
            'alignment_options': {
                'genome': genome.id,
                'mode': '--end-to-end',  # as default
                'N': 1,
                'L': 9,
                'gbar': 1,
                'mp': '2',
                'rdg': '0,2',
                'rfg': '0,2',
                'score_min': 'C,-20,0'
            },
            'quant_options': {
                'readlengths': 26,
                'alignscores': -6
            }
        }

        self.run_process('workflow-trim-align-quant', input_workflow)
        workflow = Data.objects.last()

        self.assertFile(workflow, 'rc', pf_out + 'SM31_ss_trimmed_trimmed_count_matrix.txt.gz', compression='gzip')
        self.assertFile(workflow, 'exp', pf_out + 'SM31_ss_trimmed_trimmed_count_matrix.txt.gz', compression='gzip')
        self.assertFields(workflow, 'exp_type', 'RC')
        self.assertJSON(workflow, workflow.output['exp_json'], '', pf_out + 'SM31_ss_json.txt.gz')
        self.assertFields(workflow, 'source', 'shRNA-gene-sequences')
        self.assertFields(workflow, 'species', species)
        self.assertFields(workflow, 'build', build)
        self.assertFields(workflow, 'feature_type', 'shRNA')
        self.assertFile(workflow, 'mapped_species', pf_out + 'SM31_ss_trimmed_trimmed_mapped_species.txt.gz',
                        compression='gzip')
