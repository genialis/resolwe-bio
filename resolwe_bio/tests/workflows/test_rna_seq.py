# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class RNASeqWorkflowTestCase(BioProcessTestCase):
    def test_heatseq_workflow(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz', 'source': 'DICTYBASE'}
        annotation = self.run_process('upload-gff3', inputs)

        self.run_process(
            'workflow-rnaseq-cuffquant', {
                'reads': reads.id,
                'genome': genome.id,
                'annotation': annotation.id,
            }
        )
