from .base import BaseProcessorTestCase
from .utils import PreparedData
from server.models import Data


class ClusteringProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_hc_clustering(self):
        """Cannot use assertJSON - JSON output contains ETC object IDs."""
        genome = self.prepare_genome()
        reads1 = self.prepare_reads('00Hr.fastq.gz')
        reads2 = self.prepare_reads('20Hr.fastq.gz')

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs, Data.STATUS_DONE)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads1.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_1 = self.run_processor('alignment:tophat-2-0-13', inputs, Data.STATUS_DONE)

        inputs = {
            'genome': genome.pk,
            'reads': reads2.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_2 = self.run_processor('alignment:tophat-2-0-13', inputs, Data.STATUS_DONE)

        inputs = {
            'genome': genome.pk,
            'gff': annotation.pk}
        mappability = self.run_processor('mappability:bcm-1-0-0', inputs, Data.STATUS_DONE)

        inputs = {
            'alignment': aligned_reads_1.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_1 = self.run_processor('expression:bcm-1-0-0', inputs, Data.STATUS_DONE)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_2 = self.run_processor('expression:bcm-1-0-0', inputs, Data.STATUS_DONE)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_3 = self.run_processor('expression:bcm-1-0-0', inputs, Data.STATUS_DONE)

        inputs = {'expressions': [expression_1.pk, expression_2.pk, expression_3.pk]}
        etc = self.run_processor('etc:bcm-1-0-0', inputs, Data.STATUS_DONE)

        inputs = {
            'etcs': [etc.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067100']}
        self.run_processor('clustering:hierarchical:bcm-1-0-0', inputs, Data.STATUS_DONE)
