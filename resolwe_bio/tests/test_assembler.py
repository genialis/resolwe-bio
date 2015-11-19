# pylint: disable=missing-docstring
from .utils import ProcessTestCase


class AbyssProcessorTestCase(ProcessTestCase):
    def test_abyss(self):
        se_reads = self.prepare_reads('reads.fastq.gz')

        inputs = {'src1': 'reads_paired_abyss_1.fastq.gz', 'src2': 'reads_paired_abyss_2.fastq.gz'}
        reads = self.run_processor('import:upload:reads-fastq-paired-end', inputs)

        inputs = {'reads': reads.pk,
                  'se': se_reads.pk,
                  'options': {'k': 25,
                              'name': 'ecoli'}}

        self.run_processor('assembler:abyss', inputs)
        # self.assertFiles(abyss, 'contigs', 'contigs.zip', compression='zip')
