import os
import unittest

from django.conf import settings

from resolwe.flow.models import Process, iterate_schema
from resolwe.flow.utils.test import ProcessTestCase


TEST_FILES_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'tests',
    'processes',
    'files',
)


def skipDockerFailure(reason):
    """Skip the decorated test unless ``TESTS_SKIP_DOCKER_FAILURES``
    Django setting is set to ``False``. ``reason`` should describe why
    the test is being skipped.
    """
    if getattr(settings, 'TESTS_SKIP_DOCKER_FAILURES', True):
        return unittest.skip(reason)
    return lambda func: func


def skipUnlessLargeFiles(*files):
    """Skip the decoreted tests unless the given large files are available.

    :param list *files: variable lenght files list, where each element
        represents a large file path relative to the ``TEST_FILES_DIR``
        directory

    """
    for file_ in files:
        if not os.path.isfile(os.path.join(TEST_FILES_DIR, file_)):
            return unittest.skip("File '{}' is not available".format(file_))
        try:
            with open(os.path.join(TEST_FILES_DIR, file_)) as f:
                if f.readline().startswith('version https://git-lfs.github.com/spec/'):
                    return unittest.skip("Only Git LFS pointer is available "
                                         "for file '{}'".format(file_))
        except UnicodeDecodeError:
            # file_ is a binary file (this is expected)
            pass
    return lambda func: func


class BioProcessTestCase(ProcessTestCase):
    def setUp(self):
        super(BioProcessTestCase, self).setUp()
        self.files_path = TEST_FILES_DIR

    def prepare_genome(self):
        """Prepare genome FASTA."""
        inputs = {"src": "genome.fasta.gz",
                  "bowtie_index": "bt_index.tar.gz",
                  "bowtie2_index": "bt2_index.tar.gz",
                  "bwa_index": "bwa_index.tar.gz",
                  "hisat2_index": "hisat2_index.tar.gz"}
        return self.run_processor('upload-genome', inputs)

    def prepare_reads(self, fn=['reads.fastq.gz']):
        """Prepare NGS reads FASTQ."""
        inputs = {'src': fn}
        return self.run_processor('upload-fastq-single', inputs)

    def prepare_paired_reads(self, fw=['fw_reads.fastq.gz'], rw=['rw_reads.fastq.gz']):
        """Prepare NGS reads FASTQ."""
        inputs = {'src1': fw, 'src2': rw}
        return self.run_processor('upload-fastq-paired', inputs)

    def prepare_bam(self, fn='sp_test.bam'):
        """Prepare alignment BAM."""
        inputs = {'src': fn}
        return self.run_processor('upload-bam', inputs)

    def prepare_annotation(self, fn='sp_test.gtf'):
        """Prepare annotation GTF."""
        inputs = {'src': fn}
        return self.run_processor('upload-gtf', inputs)

    def prepare_expression(self, f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM", name='Expression'):
        """Prepare expression."""
        inputs = {'rc': f_rc, 'exp': f_exp, 'exp_type': f_type, 'exp_name': name}
        return self.run_processor('upload-expression', inputs)
