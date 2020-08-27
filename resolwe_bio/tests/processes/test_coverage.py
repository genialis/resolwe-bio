from os.path import join

from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class CoverageProcessorTestCase(BioProcessTestCase):
    @skipUnlessLargeFiles("56GSID_10k_mate1_RG.bam")
    @tag_process("coveragebed")
    def test_amplicon_coverage(self):
        with self.preparation_stage():
            bam_input = {
                "src": join("large", "56GSID_10k_mate1_RG.bam"),
                "species": "Homo sapiens",
                "build": "b37",
            }
            bam = self.run_process("upload-bam", bam_input)
            master_file = self.prepare_amplicon_master_file()

        coverage = self.run_process(
            "coveragebed",
            {
                "alignment": bam.id,
                "master_file": master_file.id,
            },
        )
        self.assertFile(coverage, "cov_metrics", "56GSID_10k_covMetrics.txt")
        self.assertFile(coverage, "mean_cov", "56GSID_10k_ampmeancov.covd")
        self.assertFileExists(coverage, "amplicon_cov")
        self.assertFileExists(coverage, "covplot_html")
