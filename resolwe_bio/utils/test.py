"""Test helper functions."""

import os
import unittest

from django.conf import settings
from django.core.management import call_command
from django.test import LiveServerTestCase

from resolwe.flow.models import Collection
from resolwe.flow.models.annotations import (
    AnnotationField,
    AnnotationGroup,
    AnnotationType,
)
from resolwe.test import ProcessTestCase

TEST_FILES_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "tests", "files")
)
TEST_LARGE_FILES_DIR = os.path.join(TEST_FILES_DIR, "large")


def skipDockerFailure(reason):
    """Skip decorated tests due to failures when run in Docker.

    Unless ``TESTS_SKIP_DOCKER_FAILURES`` Django setting is set to
    ``False``. ``reason`` should describe why the test is being skipped.
    """
    if getattr(settings, "TESTS_SKIP_DOCKER_FAILURES", True):
        return unittest.skip(reason)
    return lambda func: func


def skipUnlessLargeFiles(*files):
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
                if f.readline().startswith("version https://git-lfs.github.com/spec/"):
                    return unittest.skip(
                        "Only Git LFS pointer is available "
                        "for file '{}'".format(file_path)
                    )
        except UnicodeDecodeError:
            # file_ is a binary file (this is expected)
            pass
    return lambda func: func


class BioProcessTestCase(ProcessTestCase):
    """Base class for writing bioinformatics process tests.

    It is a subclass of Resolwe's :class:`~resolwe.test.ProcessTestCase`
    with some specific functions used for testing bioinformatics
    processes.

    """

    def setUp(self):
        """Initialize test files path and species annotation."""
        super(BioProcessTestCase, self).setUp()
        self.files_path = TEST_FILES_DIR

        general_group = AnnotationGroup.objects.create(name="general", sort_order=1)
        AnnotationField.objects.create(
            name="species",
            sort_order=1,
            group=general_group,
            type=AnnotationType.STRING.value,
            vocabulary={
                "Caenorhabditis elegans": "Caenorhabditis elegans",
                "Cricetulus griseus": "Cricetulus griseus",
                "Dictyostelium discoideum": "Dictyostelium discoideum",
                "Dictyostelium purpureum": "Dictyostelium purpureum",
                "Drosophila melanogaster": "Drosophila melanogaster",
                "Homo sapiens": "Homo sapiens",
                "Macaca mulatta": "Macaca mulatta",
                "Mus musculus": "Mus musculus",
                "Rattus norvegicus": "Rattus norvegicus",
                "other": "Other",
            },
        )

        biospecimen_group = AnnotationGroup.objects.create(
            name="biospecimen_information", sort_order=2
        )
        AnnotationField.objects.create(
            name="organ",
            sort_order=1,
            group=biospecimen_group,
            type=AnnotationType.STRING.value,
        )

    def prepare_reads(self, fn=["reads.fastq.gz"]):
        """Prepare NGS reads FASTQ."""
        inputs = {"src": fn}
        return self.run_process("upload-fastq-single", inputs)

    def prepare_paired_reads(
        self, mate1=["fw reads.fastq.gz"], mate2=["rw reads.fastq.gz"]
    ):
        """Prepare NGS reads FASTQ."""
        inputs = {"src1": mate1, "src2": mate2}
        return self.run_process("upload-fastq-paired", inputs)

    def prepare_bam(
        self, fn="sp_test.bam", species="Dictyostelium discoideum", build="dd-05-2009"
    ):
        """Prepare alignment BAM."""
        inputs = {"src": fn, "species": species, "build": build}
        return self.run_process("upload-bam", inputs)

    def prepare_annotation(
        self,
        fn="sp_test.gtf",
        source="DICTYBASE",
        species="Dictyostelium discoideum",
        build="dd-05-2009",
    ):
        """Prepare annotation GTF."""
        inputs = {"src": fn, "source": source, "species": species, "build": build}
        return self.run_process("upload-gtf", inputs)

    def prepare_annotation_gff(
        self,
        fn="annotation dicty.gff.gz",
        source="DICTYBASE",
        species="Dictyostelium discoideum",
        build="dd-05-2009",
    ):
        """Prepare annotation GFF3."""
        inputs = {"src": fn, "source": source, "species": species, "build": build}
        return self.run_process("upload-gff3", inputs)

    def prepare_ref_seq(
        self, fn="adapters.fasta", species="Other", build="Illumina adapters"
    ):
        """Prepare reference sequence FASTA."""
        return self.run_process(
            "upload-fasta-nucl",
            {
                "src": fn,
                "species": species,
                "build": build,
            },
        )

    def prepare_expression(
        self,
        f_rc="exp_1_rc.tab.gz",
        f_exp="exp_1_tpm.tab.gz",
        f_type="TPM",
        name="Expression",
        source="DICTYBASE",
        feature_type="gene",
        species="Dictyostelium discoideum",
        build="dd-05-2009",
    ):
        """Prepare expression."""
        inputs = {
            "exp": f_exp,
            "exp_type": f_type,
            "exp_name": name,
            "source": source,
            "species": species,
            "build": build,
            "feature_type": feature_type,
        }
        if f_rc:
            inputs["rc"] = f_rc
        expression = self.run_process("upload-expression", inputs)
        return expression


class KBBioProcessTestCase(BioProcessTestCase, LiveServerTestCase):
    """Class for bioinformatics process tests that use knowledge base.

    It is based on :class:`~.BioProcessTestCase` and Django's
    :class:`~django.test.LiveServerTestCase`.
    The latter launches a live Django server in a separate thread so
    that the tests may use it to query the knowledge base.

    """

    # This is work-around since Django's LiveServerTestCase apparently doesn't
    # properly honor Django settings set in tests/settings.py.
    _overridden_settings = {
        # If a process that uses the live Django server fails, we want to see
        # full debugging output rather than just laconic message like
        # "Server Error (500)".
        "DEBUG": True,
    }

    def setUp(self):
        """Set up test gene information knowledge base, create collection."""
        super().setUp()

        self.collection = Collection.objects.create(
            name="Test collection", contributor=self.admin
        )

        call_command(
            "insert_features", os.path.join(TEST_FILES_DIR, "features_gsea.tab.zip")
        )
        call_command(
            "insert_mappings", os.path.join(TEST_FILES_DIR, "mappings_gsea.tab.zip")
        )

    def run_process(self, *args, **kwargs):
        """Run processes in collection."""
        kwargs["collection"] = kwargs.get("collection", self.collection)
        return super().run_process(*args, **kwargs)
