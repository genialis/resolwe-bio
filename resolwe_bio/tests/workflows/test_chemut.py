from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.filter import filter_vcf_variable
from resolwe_bio.utils.test import BioProcessTestCase


class CheMutWorkflowTestCase(BioProcessTestCase):
    @tag_process("workflow-chemut")
    def test_chemut_workflow(self):
        with self.preparation_stage():
            inputs = {
                "src": "chemut_genome.fasta.gz",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            ref_seq = self.run_process("upload-fasta-nucl", inputs)
            bwa_index = self.run_process("bwa-index", {"ref_seq": ref_seq.id})

            inputs = {"src1": ["AX4_mate1.fq.gz"], "src2": ["AX4_mate2.fq.gz"]}

            parental_reads = self.run_process("upload-fastq-paired", inputs)

            inputs = {"src1": ["CM_mate1.fq.gz"], "src2": ["CM_mate2.fq.gz"]}

            mut_reads = self.run_process("upload-fastq-paired", inputs)

            inputs = {"genome": bwa_index.id, "reads": parental_reads.id}
            align_parental = self.run_process("alignment-bwa-mem", inputs)

            inputs = {"genome": bwa_index.id, "reads": mut_reads.id}
            align_mut = self.run_process("alignment-bwa-mem", inputs)

        self.run_process(
            "workflow-chemut",
            {
                "analysis_type": "snv",
                "parental_strains": [align_parental.id],
                "mutant_strains": [align_mut.id],
                "genome": ref_seq.id,
                "Vc": {"stand_emit_conf": 15, "stand_call_conf": 35, "rf": True},
                "Vf": {"read_depth": 7},
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        variants = Data.objects.last()
        self.assertFile(
            variants,
            "vcf",
            "chemut.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )
