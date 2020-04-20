import os

from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class WgsWorkflowTestCase(BioProcessTestCase):
    @tag_process("workflow-wgs-paired")
    def test_wgs_workflow(self):
        def filter_gatkcmd(line):
            return line.startswith(b"##GATKCommandLine")

        with self.preparation_stage():
            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "./bqsr/input/hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )
            bwa_index = self.run_process("bwa-index", {"ref_seq": ref_seq.id})

            reads = self.prepare_paired_reads(
                mate1=["./workflow_wes/input/TP53_1.fastq.gz"],
                mate2=["./workflow_wes/input/TP53_2.fastq.gz"],
            )

            kbase = []
            for i in ["./bqsr/input/dbsnp_TP53.vcf.gz"]:
                kbase.append(
                    self.run_process(
                        "upload-variants-vcf",
                        {"src": i, "species": "Homo sapiens", "build": "custom_build"},
                    )
                )
            adapters = self.prepare_ref_seq()

        self.run_process(
            "workflow-wgs-paired",
            {
                "bwa_index": bwa_index.id,
                "ref_seq": ref_seq.id,
                "reads": reads.id,
                "known_sites": [i.id for i in kbase],
                "hc_dbsnp": kbase[0].id,
                "advanced": {
                    "trimming": {
                        "adapters": adapters.id,
                        "seed_mismatches": 2,
                        "simple_clip_threshold": 10,
                        "min_adapter_length": 8,
                        "palindrome_clip_threshold": 30,
                        "leading": 20,
                        "trailing": 3,
                        "minlen": 40,
                    },
                    "align": {"m": True, "scoring": {"unpaired_p": 17},},
                    "bqsr": {
                        "read_group": "-LB=DAB;-PL=Illumina;-PU=barcode;-SM=sample1"
                    },
                    "hc": {"stand_call_conf": 3, "mbq": 3},
                },
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        variants = Data.objects.filter(process__slug="vc-gatk4-hc").last()
        self.assertFile(
            variants,
            "vcf",
            os.path.join("wgs_workflow", "output", "tp53_1fastqgz.gatkHC.vcf.gz"),
            compression="gzip",
            file_filter=filter_gatkcmd,
        )
