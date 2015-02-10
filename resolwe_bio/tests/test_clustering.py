from .base import BaseProcessorTestCase


class ClusteringProcessorTestCase(BaseProcessorTestCase):
    def prepair_genome(self):
        inputs = {'src': 'genome.fasta.gz'}
        genome = self.run_processor('import:upload:genome-fasta', inputs)
        self.assertDone(genome)
        self.assertFiles(genome, 'fasta', 'genome.fasta.gz')
        return genome

    def test_hc_clustering(self):
        """Cannot use assertJSON - JSON output contains ETC object IDs."""
        genome = self.prepair_genome()

        inputs = {'src': '00Hr.fastq.gz'}
        reads1 = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(reads1)

        inputs = {'src': '20Hr.fastq.gz'}
        reads2 = self.run_processor('import:upload:reads-fastq', inputs)
        self.assertDone(reads2)

        inputs = {'src': 'annotation.gff'}
        annotation = self.run_processor('import:upload:annotation-gff3', inputs)
        self.assertDone(annotation)
        self.assertFiles(annotation, 'gff', 'annotation.gff')

        inputs = {
            'genome': genome.pk,
            'reads': reads1.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_1 = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads_1)

        inputs = {
            'genome': genome.pk,
            'reads': reads2.pk,
            'gff': annotation.pk,
            'PE_options': {
                'library_type': "fr-unstranded"}}
        aligned_reads_2 = self.run_processor('alignment:tophat-2-0-13', inputs)
        self.assertDone(aligned_reads_2)

        inputs = {
            'genome': genome.pk,
            'gff': annotation.pk}
        mappability = self.run_processor('mappability:bcm-1-0-0', inputs)
        self.assertDone(mappability)

        inputs = {
            'alignment': aligned_reads_1.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_1 = self.run_processor('expression:bcm-1-0-0', inputs)
        self.assertDone(expression_1)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_2 = self.run_processor('expression:bcm-1-0-0', inputs)
        self.assertDone(expression_2)

        inputs = {
            'alignment': aligned_reads_2.pk,
            'gff': annotation.pk,
            'mappable': mappability.pk}
        expression_3 = self.run_processor('expression:bcm-1-0-0', inputs)
        self.assertDone(expression_3)

        inputs = {'expressions': [expression_1.pk, expression_2.pk, expression_3.pk]}
        etc = self.run_processor('etc:bcm-1-0-0', inputs)
        self.assertDone(etc)

        inputs = {
            'etcs': [etc.pk],
            'genes': ['DPU_G0067096', 'DPU_G0067098', 'DPU_G0067100']}
        clustering = self.run_processor('clustering:hierarchical:bcm-1-0-0', inputs)
        self.assertDone(clustering)
