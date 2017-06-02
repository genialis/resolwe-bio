# pylint: disable=missing-docstring
from os.path import join

from resolwe_bio.utils.test import skipUnlessLargeFiles, BioProcessTestCase


class ReportProcessorTestCase(BioProcessTestCase):

    @skipUnlessLargeFiles('56GSID_10k_mate1_RG.bam')
    def test_amplicon_report(self):
        template = self.run_process('upload-file', {'src': 'report_template.tex'})
        logo = self.run_process('upload-file', {'src': 'genialis_logo.pdf'})

        bam = self.run_process('upload-bam', {'src': join('large', '56GSID_10k_mate1_RG.bam')})
        master_file = self.run_process('upload-master-file', {'src': '56G_masterfile_test.txt'})
        coverage = self.run_process('coveragebed', {'alignment': bam.id, 'master_file': master_file.id})

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

        report_inputs = {
            'coverage': coverage.id,
            'pcr_metrics': pcr_metrics.id,
            'template': template.id,
            'logo': logo.id,
            'panel_name': '56G Oncology Panel v2',
            'annot_vars': [annot_variants.id]
        }

        self.run_process('amplicon-report', report_inputs)
