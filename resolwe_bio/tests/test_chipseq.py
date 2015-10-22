# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class ChipSeqProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_chipseq(self):
        inputs = {"src": "ChIP-Seq-Control.bam"}
        control_bam = self.run_processor("import:upload:mapping-bam", inputs)

        inputs = {"src": "ChIP-Seq-Case.bam"}
        case_bam = self.run_processor("import:upload:mapping-bam", inputs)

        inputs = {"src": "dd_genes_chr1.bed"}
        bed = self.run_processor("import:upload:bed", inputs)

        inputs = {
            'case': case_bam.pk,
            'control': control_bam.pk,
            'settings': {'nomodel': True,
                         'gsize': '3.4e7',
                         'pvalue': 0.001,
                         'slocal': 2000,
                         'extsize': 100,
                         'call_summits': True}}
        macs2 = self.run_processor("macs2:callpeak", inputs)

        inputs = {
            'peaks': macs2.pk,
            'bed': bed.pk}
        peak_score = self.run_processor("chipseq:peakscore", inputs)
        self.assertFiles(peak_score, "peak_score", "macs2_peakscore_genomicContext")

        inputs = {"peakscore": peak_score.pk}
        gene_score = self.run_processor("chipseq:genescore", inputs)
        self.assertFiles(gene_score, "genescore", "macs2_geneScore.xls")
