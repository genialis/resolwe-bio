# pylint: disable=missing-docstring
from .utils import skipDockerFailure, BioProcessTestCase


class PseudoalignmentProcessorTestCase(BioProcessTestCase):

    @skipDockerFailure("Fails with: kallisto: command not found")
    def test_kallisto(self):
        inputs = {"src": "cDNA_reference.fasta.gz"}
        cdna = self.run_processor("upload-fasta-nucl", inputs)

        inputs = {"src": "reads_pseudoaligner.fastq.gz"}
        reads = self.run_processor("upload-fastq-single", inputs)

        inputs = {"fasta": cdna.pk}
        kallisto_index = self.run_processor("kallisto-index", inputs)

        inputs = {
            'reads': reads.pk,
            'transcriptome': kallisto_index.pk,
            'options': {
                'single': True,
                'fragment_length_single': 75}}
        kallisto = self.run_processor("kallisto-quant", inputs)
        self.assertFiles(kallisto, "abundance", "kallisto_abundance.txt.gz", compression='gzip')
        self.assertFiles(kallisto, "rc", "kallisto_abundance_rc.tab.gz", compression='gzip')
