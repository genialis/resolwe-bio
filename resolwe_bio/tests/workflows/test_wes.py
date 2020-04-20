from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class WESTestCase(BioProcessTestCase):
    @tag_process("workflow-wes")
    def test_wes(self):
        def filter_gatkcmd(line):
            return line.startswith(b"##GATKCommandLine")

        species = "Homo sapiens"
        build = "custom_build"

        with self.preparation_stage():
            reads = self.prepare_paired_reads(
                mate1=["./workflow_wes/input/TP53_1.fastq.gz"],
                mate2=["./workflow_wes/input/TP53_2.fastq.gz"],
            )

            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "./bqsr/input/hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": species,
                    "build": build,
                },
            )

            bwa_index = self.run_process("bwa-index", {"ref_seq": ref_seq.id})

            bc_bedpe = self.run_process(
                "upload-bedpe",
                {
                    "src": "./bamclipper/input/TP53.bedpe",
                    "species": species,
                    "build": build,
                },
            )

            kbase = []
            for i in ["./bqsr/input/dbsnp_TP53.vcf.gz"]:
                kbase.append(
                    self.run_process(
                        "upload-variants-vcf",
                        {"src": i, "species": species, "build": build},
                    )
                )

            intervals = self.run_process(
                "upload-bed",
                {"src": "./bqsr/input/TP53.bed", "species": species, "build": build},
            )

            adapters = self.prepare_ref_seq()

        input_workflow = {
            "reads": reads.pk,
            "bwa_index": bwa_index.id,
            "ref_seq": ref_seq.id,
            "known_sites": [i.id for i in kbase],
            "intervals": intervals.id,
            "hc_dbsnp": kbase[0].id,
            "validation_stringency": "LENIENT",
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
                "align": {
                    "seed_l": 19,
                    "band_w": 100,
                    "m": True,
                    "re_seeding": 1.5,
                    "scoring": {
                        "match": 1,
                        "mismatch": 4,
                        "gap_o": 6,
                        "gap_e": 1,
                        "clipping": 5,
                        "unpaired_p": 17,
                    },
                    "report_tr": 30,
                },
                "bamclipper": {"bedpe": bc_bedpe.id, "skip": False},
                "markduplicates": {
                    "md_skip": False,
                    "md_remove_duplicates": False,
                    "md_assume_sort_order": "",
                },
                "bqsr": {"read_group": "-LB=DAB;-PL=Illumina;-PU=barcode;-SM=sample1"},
                "hc": {"stand_call_conf": 3, "mbq": 3},
            },
        }

        self.run_process("workflow-wes", input_workflow)
        wes = Data.objects.last()
        self.assertFile(
            wes,
            "vcf",
            "./workflow_wes/output/tp53_1fastqgz.gatkHC.vcf.gz",
            compression="gzip",
            file_filter=filter_gatkcmd,
        )

        # Test skipping bamclipper
        input_workflow_skipbc = input_workflow
        input_workflow_skipbc["advanced"]["bamclipper"]["skip"] = True
        self.run_process("workflow-wes", input_workflow_skipbc)
        wes_skip = Data.objects.filter(process__slug="bamclipper").last()

        self.assertEqual(wes_skip.process_info, ["Skipping bamclipper step."])
