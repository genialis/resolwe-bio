# pylint: disable=no-member
from server.models import Data


class PreparedData(object):
    """MixIn with functions for processing common data types.

    This class is MixIn meant for use with
    :class:`BaseProcessorTestCase` and includes functions for processing
    common data types as inputs for more advanced tests.

    """
    def prepare_genome(self, fn='genome.fasta.gz'):
        inputs = {'src': fn}
        return self.run_processor('import:upload:genome-fasta', inputs, Data.STATUS_DONE)

    def prepare_reads(self, fn='reads.fastq.gz'):
        inputs = {'src': fn}
        return self.run_processor('import:upload:reads-fastq', inputs, Data.STATUS_DONE)

    def prepare_bam(self, fn='sp_test.bam'):
        inputs = {'src': fn}
        return self.run_processor('import:upload:mapping-bam', inputs, Data.STATUS_DONE)

    def prepare_annotation(self, fn='sp_test.gtf'):
        inputs = {'src': fn}
        return self.run_processor('import:upload:annotation-gtf', inputs, Data.STATUS_DONE)
