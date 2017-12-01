# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase
from resolwe.flow.models import Data
from resolwe.test import tag_process


class AmpliconWorkflowTestCase(BioProcessTestCase):
    @skipDockerFailure("Processor requires a custom Docker image.")
    @tag_process('workflow-accel')
    def test_amplicon_workflow(self):
        with self.preparation_stage():
            inputs = {
                'src1': ['56GSID_10k_mate1.fastq.gz'],
                'src2': ['56GSID_10k_mate2.fastq.gz']}
            reads = self.run_process('upload-fastq-paired', inputs)

            adapters = self.run_process('upload-fasta-nucl', {'src': 'adapters.fasta'})
            genome = self.run_process('upload-genome', {'src': 'hs_b37_chr2_small.fasta.gz'})

            master_file = self.prepare_amplicon_master_file()

            inputs = {'src': '1000G_phase1.indels.b37_chr2_small.vcf.gz'}
            indels = self.run_process('upload-variants-vcf', inputs)

            dbsnp = self.run_process('upload-variants-vcf', {'src': 'dbsnp_138.b37.chr2_small.vcf.gz'})

            template = self.run_process('upload-file', {'src': 'report_template.tex'})
            logo = self.run_process('upload-file', {'src': 'genialis_logo.pdf'})
            template_html = self.run_process('upload-file', {'src': 'report_html_template.html'})
            bokeh_css = self.run_process('upload-file', {'src': 'bokeh-0.12.9.min.css'})
            bokeh_js = self.run_process('upload-file', {'src': 'bokeh-0.12.9.min.js'})

        self.run_process(
            'workflow-accel', {
                'reads': reads.id,
                'genome': genome.id,
                'master_file': master_file.id,
                'adapters': adapters.id,
                'template_html': template_html.id,
                'bokeh_css': bokeh_css.id,
                'bokeh_js': bokeh_js.id,
                'preprocess_bam': {
                    'known_vars': [dbsnp.id],
                    'known_indels': [indels.id]
                },
                'gatk': {
                    'dbsnp': dbsnp.id,
                    'mbq': 20,
                    'stand_emit_conf': 20
                },
                'lofreq': {
                    'min_bq': 20,
                    'min_alt_bq': 20
                },
                'var_annot': {
                    'known_vars_db': [dbsnp.id]
                },
                'report': {
                    'template': template.id,
                    'logo': logo.id,
                },
                'threads': 2
            }
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
