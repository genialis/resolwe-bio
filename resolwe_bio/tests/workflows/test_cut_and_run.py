from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class CutAndRunTestCase(BioProcessTestCase):
    @tag_process("workflow-cutnrun")
    def test_cutnrun_workflow(self):
        species = "Homo sapiens"
        build = "GRCh38_ens109"

        with self.preparation_stage():
            # Data is from chr1:1-1000000 of human sample.
            reads = self.prepare_paired_reads(
                mate1=["./workflow_cutnrun/input/mate1.fq.gz"],
                mate2=["./workflow_cutnrun/input/mate2.fq.gz"],
            )

            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    # Human genome prepared by cutting chr1 up to 1M bp.
                    "src": "./workflow_cutnrun/input/GRCh38_chr1_1_1M.fasta.gz",
                    "species": species,
                    "build": build,
                },
            )

            bowtie2_index = self.run_process("bowtie2-index", {"ref_seq": ref_seq.id})

            promoter_regions = self.run_process(
                "upload-bed",
                {
                    "src": "./workflow_cutnrun/input/promoters_hsapiens_grch38_ensembl_109_small.bed",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens109",
                },
            )

        input_workflow = {
            "reads": reads.id,
            "peak_calling_options": {
                "promoter": promoter_regions.id,
            },
            "alignment_options": {
                "genome": bowtie2_index.id,
                "spikein_genome": bowtie2_index.id,
            },
        }

        self.run_process("workflow-cutnrun", input_workflow)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        cnr = Data.objects.filter(process__slug="calculate-bigwig").last()
        self.assertFile(
            obj=cnr,
            field_path="bigwig",
            fn="./workflow_cutnrun/output/normalized.bigwig",
        )

        multiqc = Data.objects.filter(process__slug="multiqc").last()
        self.assertFileExists(multiqc, "report")
