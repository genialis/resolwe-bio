# pylint: disable=missing-docstring
from resolwe_bio.utils.test import skipDockerFailure, BioProcessTestCase


class DiffExpProcessorTestCase(BioProcessTestCase):

    def test_cuffdiff(self):
        inputs = {"src": "cuffquant_1.cxb"}
        cuffquant = self.run_processor("upload-cxb", inputs)

        inputs = {"src": "cuffquant_2.cxb"}
        cuffquant2 = self.run_processor("upload-cxb", inputs)

        annotation = self.prepare_annotation(fn='hg19_chr20_small.gtf.gz')

        inputs = {
            'case': [cuffquant.id],
            'control': [cuffquant2.id],
            'annotation': annotation.id}
        cuffdiff = self.run_processor('cuffdiff', inputs)
        self.assertFile(cuffdiff, 'diffexp', 'cuffdiff.tab.gz', compression='gzip')
        self.assertJSON(cuffdiff, cuffdiff.output['de_data'], '', 'cuffdiff.json.gz')

    def test_bayseq_bcm(self):
        expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")

        mappa = self.run_processor("upload-mappability", {"src": "purpureum_mappability_50.tab.gz"})

        inputs = {
            'name': "00vs20",
            'case': [expression_1.id],
            'control': [expression_2.id],
            'replicates': ['1', '2'],
            'mappability': mappa.id}
        diff_exp = self.run_processor('differentialexpression-bcm', inputs)
        self.assertJSON(diff_exp, diff_exp.output['volcano_plot'], '', 'bayseq_volcano.json.gz')

    def test_deseq2(self):
        expression_1 = self.prepare_expression(f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM")
        expression_2 = self.prepare_expression(f_rc='exp_2_rc.tab.gz', f_exp='exp_2_tpm.tab.gz', f_type="TPM")

        inputs = {
            'case': [expression_1.pk],
            'control': [expression_2.pk]
        }

        diff_exp = self.run_processor('differentialexpression-deseq2', inputs)
        self.assertFile(diff_exp, "diffexp", 'diffexp_deseq2.tab.gz', compression='gzip')

    def test_limma(self):
        expression_1 = self.prepare_expression(f_exp='exp_limma_1.tab.gz', f_type="Log2")
        expression_2 = self.prepare_expression(f_exp='exp_limma_2.tab.gz', f_type="Log2")
        expression_3 = self.prepare_expression(f_exp='exp_limma_3.tab.gz', f_type="Log2")
        expression_4 = self.prepare_expression(f_exp='exp_limma_4.tab.gz', f_type="Log2")

        inputs = {
            'case': [expression_1.pk, expression_2.pk],
            'control': [expression_3.pk, expression_4.pk]
        }

        diff_exp = self.run_processor('differentialexpression-limma', inputs)
        self.assertFile(diff_exp, "diffexp", 'diffexp_limma.tab.gz', compression='gzip')
