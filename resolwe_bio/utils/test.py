"""Test helper functions."""
import os
import unittest

from django.conf import settings

from resolwe.flow.utils.test import ProcessTestCase

TEST_FILES_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'tests', 'files'))
TEST_LARGE_FILES_DIR = os.path.join(TEST_FILES_DIR, 'large')


def skipDockerFailure(reason):  # pylint: disable=invalid-name
    """Skip decorated tests due to failures when run in Docker.

    Unless ``TESTS_SKIP_DOCKER_FAILURES`` Django setting is set to
    ``False``. ``reason`` should describe why the test is being skipped.
    """
    if getattr(settings, 'TESTS_SKIP_DOCKER_FAILURES', True):
        return unittest.skip(reason)
    return lambda func: func


def skipUnlessLargeFiles(*files):  # pylint: disable=invalid-name
    r"""Skip decorated tests unless large files are available.

    :param list \*files: variable lenght files list, where each element
                         represents a large file path relative to the
                         ``TEST_LARGE_FILES_DIR`` directory

    """
    for file_ in files:
        file_path = os.path.join(TEST_LARGE_FILES_DIR, file_)
        if not os.path.isfile(file_path):
            return unittest.skip("File '{}' is not available".format(file_path))
        try:
            with open(file_path) as f:
                if f.readline().startswith('version https://git-lfs.github.com/spec/'):
                    return unittest.skip("Only Git LFS pointer is available "
                                         "for file '{}'".format(file_path))
        except UnicodeDecodeError:
            # file_ is a binary file (this is expected)
            pass
    return lambda func: func


class BioProcessTestCase(ProcessTestCase):
    """Test class for bioinformatics processes."""

    def setUp(self):
        """Initialize test files path."""
        super(BioProcessTestCase, self).setUp()
        self.files_path = TEST_FILES_DIR

    def prepare_genome(self):
        """Prepare genome FASTA."""
        inputs = {"src": "genome.fasta.gz",
                  "bowtie_index": "bt_index.tar.gz",
                  "bowtie2_index": "bt2_index.tar.gz",
                  "bwa_index": "bwa_index.tar.gz",
                  "hisat2_index": "hisat2_index.tar.gz"}
        return self.run_process('upload-genome', inputs)

    def prepare_reads(self, fn=['reads.fastq.gz']):
        """Prepare NGS reads FASTQ."""
        inputs = {'src': fn}
        return self.run_process('upload-fastq-single', inputs)

    def prepare_paired_reads(self, mate1=['fw_reads.fastq.gz'], mate2=['rw_reads.fastq.gz']):
        """Prepare NGS reads FASTQ."""
        inputs = {'src1': mate1, 'src2': mate2}
        return self.run_process('upload-fastq-paired', inputs)

    def prepare_bam(self, fn='sp_test.bam'):
        """Prepare alignment BAM."""
        inputs = {'src': fn}
        return self.run_process('upload-bam', inputs)

    def prepare_annotation(self, fn='sp_test.gtf'):
        """Prepare annotation GTF."""
        inputs = {'src': fn}
        return self.run_process('upload-gtf', inputs)

    def prepare_expression(self, f_rc='exp_1_rc.tab.gz', f_exp='exp_1_tpm.tab.gz', f_type="TPM", name='Expression'):
        """Prepare expression."""
        inputs = {'rc': f_rc, 'exp': f_exp, 'exp_type': f_type, 'exp_name': name}
        return self.run_process('upload-expression', inputs)
