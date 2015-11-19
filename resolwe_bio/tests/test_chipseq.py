# pylint: disable=missing-docstring
from .utils import ProcessTestCase


class ChipSeqProcessorTestCase(ProcessTestCase):
    def test_chipseq(self):
        inputs = {"src": "chip_seq_control.bam"}
        control_bam = self.run_processor("import:upload:mapping-bam", inputs)

        inputs = {"src": "chip_seq_case.bam"}
        case_bam = self.run_processor("import:upload:mapping-bam", inputs)

        inputs = {"src": "chip_seq.bed"}
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
        self.assertFiles(peak_score, "peak_score", "chip_seq_peakscore_genomicContext")

        inputs = {"peakscore": peak_score.pk}
        gene_score = self.run_processor("chipseq:genescore", inputs)
        self.assertFiles(gene_score, "genescore", "chip_seq_geneScore.xls")
