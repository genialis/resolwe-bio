"""Example workflow test."""

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase

class DocsProcessTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("upload-fastq-paired-docs", "workflow-qc-docs")
    def test_qc_workflow(self):

        # place the test files in the
        # resolwe_bio/tests/files/ directory
        mate1 = "mate1.fastq.gz" # replace with actual test file
        mate2 = "mate2.fastq.gz" # replace with actual test file

        with self.preparation_stage():
            fastq_reads = self.run_process(
                process_slug="upload-fastq-paired-docs",
                input_ = {
                    "mate1": mate1,
                    "mate2": mate2,
                    "species": "Homo sapiens",
                },
            )

        workflow_inputs = {"reads": fastq_reads.id}
        self.run_process("workflow-qc-docs", input_=workflow_inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        fastqc = Data.objects.filter(process__slug="fastqc-paired-end").last()
        self.assertFile(fastqc, "fastqc_report_mate1", "fastqc_report_mate1.html")
