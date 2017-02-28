# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class ReportProcessorTestCase(BioProcessTestCase):

    @skipDockerFailure("Processor requires a custom Docker image.")
    def test_amplicon_report(self):
        template = self.run_process('upload-file', {'src': 'report_template.tex'})
        logo = self.run_process('upload-file', {'src': 'genialis_logo.pdf'})

        bam = self.run_process('upload-bam', {'src': '56GSID_10k_trimmed.bam'})
        bed = self.run_process('upload-bed', {'src': '56g_targets_small.bed'})

        coverage = self.run_process('coveragebed', {'alignment': bam.id, 'bed': bed.id})

        genome = self.run_process('upload-genome', {'src': 'hs_b37_chr22_frag.fasta.gz'})
        bed_picard = self.run_process('upload-bed', {'src': '56g_targets_picard_small.bed'})

        inputs = {'src': 'Mills_and_1000G_gold_standard.indels.b37.chr22_small.vcf.gz'}
        indels = self.run_process('upload-variants-vcf', inputs)

        dbsnp = self.run_process('upload-variants-vcf', {'src': 'dbsnp_138.b37.chr22_small.vcf.gz'})

        inputs = {
            'alignment': bam.id,
            'bed': bed_picard.id,
            'genome': genome.id,
            'known_indels': [indels.id],
            'known_vars': [dbsnp.id]
        }

        preprocess_bam = self.run_process('vc-preprocess-bam', inputs)

        report_inputs = {
            'bam': preprocess_bam.id,
            'coverage': coverage.id,
            'template': template.id,
            'logo': logo.id
        }

        self.run_process('amplicon-report', report_inputs)
