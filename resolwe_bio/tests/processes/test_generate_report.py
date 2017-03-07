# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class ReportProcessorTestCase(BioProcessTestCase):

    def test_amplicon_report(self):
        template = self.run_process('upload-file', {'src': 'report_template.tex'})
        logo = self.run_process('upload-file', {'src': 'genialis_logo.pdf'})

        inputs = {
            'src': '56GSID_10k_trimmed.bam.realigned.bqsrCal.bam',
            'target_pcr_metrics': '56GSID_10k_trimmed.bam.targetPCRmetrics.txt',
            'target_coverage': '56GSID_10k_trimmed.bam.perTargetCov.txt'
        }

        preprocess_bam = self.run_process('upload-bam-vc', inputs)

        inputs = {
            'cov': '56GSID_10k_trimmed.cov',
            'covd': '56GSID_10k_trimmed.covd'
        }

        coverage = self.run_process('upload-coverage', inputs)

        inputs = {
            'annotation': '56GSID.lf.finalvars.txt',
            'summary': '56GSID_1k.gatkHC_snpEff_summary.html',
            'snpeff_genes': '56GSID_1k.gatkHC_snpEff_genes.txt'
        }

        annot_variants = self.run_process('upload-snpeff', inputs)

        report_inputs = {
            'bam': preprocess_bam.id,
            'coverage': coverage.id,
            'template': template.id,
            'logo': logo.id,
            'panel_name': '56G Oncology Panel v2',
            'annot_vars': [annot_variants.id]
        }

        self.run_process('amplicon-report', report_inputs)
