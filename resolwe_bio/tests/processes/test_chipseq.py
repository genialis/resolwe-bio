# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class ChipSeqProcessorTestCase(BioProcessTestCase):

    def test_chipseq(self):
        inputs = {"src": "chip_seq_control.bam"}
        control_bam = self.run_processor("upload-bam", inputs)

        inputs = {"src": "chip_seq_case.bam"}
        case_bam = self.run_processor("upload-bam", inputs)

        inputs = {"src": "chip_seq.bed"}
        bed = self.run_processor("upload-bed", inputs)

        inputs = {
            'case': case_bam.pk,
            'control': control_bam.pk,
            'settings': {'nomodel': True,
                         'gsize': '3.4e7',
                         'pvalue': 0.001,
                         'slocal': 2000,
                         'extsize': 100,
                         'call_summits': True}}
        macs2 = self.run_processor("macs2-callpeak", inputs)

        inputs = {
            'peaks': macs2.pk,
            'bed': bed.pk}
        peak_score = self.run_processor("chipseq-peakscore", inputs)
        self.assertFile(peak_score, "peak_score", "chip_seq_peakscore_genomicContext")

        inputs = {"peakscore": peak_score.pk}
        gene_score = self.run_processor("chipseq-genescore", inputs)
        self.assertFile(gene_score, "genescore", "chip_seq_geneScore.xls")

    def test_macs14(self):
        inputs = {"src": "chip_seq_control.bam"}
        control_bam = self.run_processor("upload-bam", inputs)

        inputs = {"src": "chip_seq_case.bam"}
        case_bam = self.run_processor("upload-bam", inputs)

        inputs = {"t": case_bam.id,
                  "c": control_bam.id}
        macs14 = self.run_processor("macs14", inputs)

    @skipUnlessLargeFiles('large/rose2_case.bam', 'large/rose2_control.bam')
    def test_rose2(self):
        inputs = {'src': 'large/rose2_case.bam'}
        bam = self.run_processor('upload-bam', inputs)

        inputs = {"src": "large/rose2_control.bam"}
        control = self.run_processor("upload-bam", inputs)

        inputs = {'src': 'macs14_chr22.bed'}
        macs_peaks = self.run_processor('upload-bed', inputs)

        inputs = {"g": 'HG19',
                  "i_upload":  macs_peaks.id,
                  "r": bam.id,
                  "c": control.id,
                  "t": 2500}
        rose2 = self.run_processor("rose2", inputs)

        # remove changing lines from the rose2 output
        def filter_created(line):
            return line.startswith(b'#Created')

        self.assertFile(rose2, 'all_enhancers', 'rose2_enhancer_table.txt', filter=filter_created)
