# pylint: disable=missing-docstring
from .base import BaseProcessorTestCase
from .utils import PreparedData


class AbyssProcessorTestCase(BaseProcessorTestCase, PreparedData):
    def test_abyss(self):
        se_reads = self.prepare_reads('20Hr.fastq.gz')

        inputs = {'src1': 'abyss_reads1.fastq.gz', 'src2': 'abyss_reads2.fastq.gz'}
        reads = self.run_processor('import:upload:reads-fastq-paired-end', inputs)

        inputs = {'reads': reads.pk,
                  'se': se_reads.pk,
                  'options': {'k': 25,
                              'name': 'ecoli'}}

        self.run_processor('assembler:abyss', inputs)
        # self.assertFiles(abyss, 'contigs', 'contigs.zip', compression='zip')
