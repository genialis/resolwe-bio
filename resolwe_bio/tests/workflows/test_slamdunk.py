import os

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class SlamdunkWorkflowTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("workflow-slamdunk-paired")
    def test_slamdunk_workflow(self):
        with self.preparation_stage():
            paired_reads = self.prepare_paired_reads(
                ["hs_slamseq_R1_complemented.fastq.gz"], ["hs_slamseq_R2.fastq.gz"]
            )
            transcripts = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": os.path.join("slamseq", "input", "hs_transcript.fasta"),
                    "species": "Homo sapiens",
                    "build": "Gencode 32",
                },
            )
            bedfile = self.run_process(
                "upload-bed",
                {
                    "src": os.path.join("slamseq", "input", "hs_transcript.bed"),
                    "species": "Homo sapiens",
                    "build": "Gencode 32",
                },
            )

        self.run_process(
            "workflow-slamdunk-paired",
            {
                "reads": paired_reads.id,
                "ref_seq": transcripts.id,
                "regions": bedfile.id,
                "options": {"read_length": 75,},
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        collapsed_tcount = Data.objects.filter(process__slug="alleyoop-collapse")[0]
        self.assertFile(
            collapsed_tcount,
            "tcount",
            os.path.join("slamseq", "output", "hs_slamseq_tcount_collapsed.txt"),
        )

        tc_exp = Data.objects.filter(process__slug="slam-count")[0]
        self.assertFile(
            tc_exp,
            "exp",
            os.path.join("slamseq", "output", "tcount_exp_tpm.txt.gz"),
            compression="gzip",
        )
