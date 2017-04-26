# pylint: disable=missing-docstring
from resolwe.flow.models import Data
from resolwe_bio.utils.test import BioProcessTestCase


class IcountWorkflowTestCase(BioProcessTestCase):
    def test_icount_primary_analysis(self):
        inputs = {'src': 'icount.workflow.in.fasta.gz'}
        genome = self.run_process('upload-fasta-nucl', inputs)
        inputs = {'src': 'icount.workflow.in.gtf', 'source': 'ENSEMBL'}
        ann = self.run_process('upload-gtf', inputs)
        inputs = {'annotation': ann.pk, 'genome': genome.pk}
        segmentation = self.run_process('icount-segment', inputs)

        inputs = {'genome2': genome.pk, 'annotation': ann.pk}
        star_index = self.run_process('alignment-star-index', inputs)

        inputs = {'src': ['icount.workflow.in.fastq']}
        reads = self.run_process('upload-fastq-single', inputs)

        self.run_process('workflow-icount', {
            'reads': reads.id,
            'index': star_index.id,
            'segmentation': segmentation.id,
        })

        peaks = Data.objects.last()
        self.assertFile(peaks, "scores", "icount.workflow.out.tsv.gz", compression='gzip')

    def test_icount_demultiplex(self):
        # TODO add star-index and annotation objects upload to test ICount workflow functionality
        # TODO assert workflow output results
        self.run_process('workflow-icount-demultiplex', {
            'reads': ['iCount_test_file_small.fastq.gz'],
            'icount_annotation': 'icount_annotation.xlsx'
        })
