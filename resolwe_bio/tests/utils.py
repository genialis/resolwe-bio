import os

try:
    from resolwe.flow.tests import ProcessTestCase
except ImportError:
    # Backward compatibility to run tests on our old platform
    from server.tests.utils import ProcessTestCase

try:
    from resolwe.flow.models import Process, iterate_schema
except ImportError:
    from server.models import Processor as Process, iterate_schema


class BioProcessTestCase(ProcessTestCase):
    def setUp(self):
        super(BioProcessTestCase, self).setUp()
        self.files_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'files')

    def prepare_genome(self, fn='genome.fasta.gz'):
        """Prepare genome FASTA."""
        inputs = {'src': fn}
        return self.run_processor('import:upload:genome-fasta', inputs)

    def prepare_reads(self, fn='reads.fastq.gz'):
        """Prepare NGS reads FASTQ."""
        inputs = {'src': fn}
        return self.run_processor('import:upload:reads-fastq', inputs)

    def prepare_bam(self, fn='sp_test.bam'):
        """Prepare alignment BAM."""
        inputs = {'src': fn}
        return self.run_processor('import:upload:mapping-bam', inputs)

    def prepare_annotation(self, fn='sp_test.gtf'):
        """Prepare annotation GTF."""
        inputs = {'src': fn}
        return self.run_processor('import:upload:annotation-gtf', inputs)

    def prepare_expression(self, f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM"):
        """Prepare expression."""
        inputs = {'rc': f_rc, 'exp': f_exp, 'exp_type': f_type}
        return self.run_processor('import:upload:expression', inputs)
