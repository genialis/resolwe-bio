# pylint: disable=missing-docstring
from resolwe_bio.utils.test import BioProcessTestCase


class ChipSeqWorkflowTestCase(BioProcessTestCase):
    def test_chip_seq_workflow(self):
        genome = self.prepare_genome()
        reads = self.prepare_reads()

        inputs = {'src': 'annotation.gff.gz', 'source': 'DICTYBASE'}
        annotation = self.run_process('upload-gff3', inputs)

        self.run_process(
            'workflow-chip-seq', {
                'reads': reads.id,
                'genome': genome.id,
                'annotation': annotation.id,
                'macs_gsize': '1.2e8',
                'rose_genome': 'MM9'
            }
        )
