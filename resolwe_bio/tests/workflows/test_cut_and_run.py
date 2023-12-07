from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class CutAndRunTestCase(BioProcessTestCase):
    @tag_process("workflow-cutnrun")
    def test_cutnrun(self):
        species = "Homo sapiens"
        build = "custom_build"

        with self.preparation_stage():
            # Data is from chr17:3020000-3040800 of mouse genome (mm10).
            reads = self.prepare_paired_reads(
                mate1=["./workflow_cutnrun/input/chr17_cutnrun_1.fastq.gz"],
                mate2=["./workflow_cutnrun/input/chr17_cutnrun_2.fastq.gz"],
            )

            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    # Mouse genome (mm10) prepared by cutting chr17 up to 30040800 bp.
                    "src": "./workflow_cutnrun/input/mm10_chr17_upto_3040800.fasta.gz",
                    "species": species,
                    "build": build,
                },
            )

            bowtie2_index = self.run_process("bowtie2-index", {"ref_seq": ref_seq.id})

        input_workflow = {
            "reads": reads.id,
            "options_aln_species": {
                "genome": bowtie2_index.id,
            },
            "options_aln_spikein": {
                "genome": bowtie2_index.id,
            },
            "options_sieve": {
                "max_frag_length": 120,
            },
        }

        self.run_process("workflow-cutnrun", input_workflow)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        cnr = Data.objects.last()
        self.assertFile(
            obj=cnr,
            field_path="bigwig",
            fn="./workflow_cutnrun/output/normalized.bigwig",
        )

    @tag_process("workflow-cutnrun-beta")
    def test_cutnrun_beta(self):
        species = "Homo sapiens"
        build = "custom_build"

        with self.preparation_stage():
            # Data is from chr17:3020000-3040800 of mouse genome (mm10).
            reads = self.prepare_paired_reads(
                mate1=["./workflow_cutnrun/input/chr17_cutnrun_1.fastq.gz"],
                mate2=["./workflow_cutnrun/input/chr17_cutnrun_2.fastq.gz"],
            )

            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    # Mouse genome (mm10) prepared by cutting chr17 up to 30040800 bp.
                    "src": "./workflow_cutnrun/input/mm10_chr17_upto_3040800.fasta.gz",
                    "species": species,
                    "build": build,
                },
            )

            bowtie2_index = self.run_process("bowtie2-index", {"ref_seq": ref_seq.id})

        input_workflow = {
            "reads": reads.id,
            "alignment_options": {
                "genome": bowtie2_index.id,
                "spike_in_genome": bowtie2_index.id,
            },
        }

        self.run_process("workflow-cutnrun-beta", input_workflow)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        cnr = Data.objects.last()
        self.assertFile(
            obj=cnr,
            field_path="bigwig",
            fn="./workflow_cutnrun/output/normalized_2.bigwig",
        )
