# pylint: disable=missing-docstring
from os.path import join

from resolwe_bio.utils.test import skipUnlessLargeFiles, BioProcessTestCase


class ReportProcessorTestCase(BioProcessTestCase):

    @skipUnlessLargeFiles('56GSID_10k_mate1_RG.bam')
    def test_amplicon_report(self):
        bam = self.run_process('upload-bam', {'src': join('large', '56GSID_10k_mate1_RG.bam')})
        master_file = self.prepare_amplicon_master_file()
        template_html = self.run_process('upload-file', {'src': 'report_html_template.html'})
        bokeh_css = self.run_process('upload-file', {'src': 'bokeh-0.12.9.min.css'})
        bokeh_js = self.run_process('upload-file', {'src': 'bokeh-0.12.9.min.js'})

        coverage = self.run_process('coveragebed', {
            'alignment': bam.id,
            'master_file': master_file.id,
            'template_html': template_html.id,
            'bokeh_css': bokeh_css.id,
            'bokeh_js': bokeh_js.id,
        })

        inputs = {
            'target_pcr_metrics': '56gsid_10k.targetPCRmetrics.txt',
            'target_coverage': '56gsid_10k.perTargetCov.txt'
        }
        pcr_metrics = self.run_process('upload-picard-pcrmetrics', inputs)

        inputs = {
            'annotation': '56GSID.lf.finalvars.txt',
            'summary': '56GSID_1k.gatkHC_snpEff_summary.html',
            'snpeff_genes': '56GSID_1k.gatkHC_snpEff_genes.txt'
        }
        annot_variants = self.run_process('upload-snpeff', inputs)

        template = self.run_process('upload-file', {'src': 'report_template.tex'})
        logo = self.run_process('upload-file', {'src': 'genialis_logo.pdf'})
        report_inputs = {
            'coverage': coverage.id,
            'pcr_metrics': pcr_metrics.id,
            'template': template.id,
            'logo': logo.id,
            'master_file': master_file.id,
            'annot_vars': [annot_variants.id]
        }
        self.run_process('amplicon-report', report_inputs)
