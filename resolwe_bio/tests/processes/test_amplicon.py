from os.path import join

from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class AmpliconProcessorTestCase(BioProcessTestCase):
    @tag_process("align-bwa-trim")
    def test_bwa_trim(self):
        with self.preparation_stage():
            inputs = {
                "src1": ["56GSID_10k_mate1.fastq.gz"],
                "src2": ["56GSID_10k_mate2.fastq.gz"],
            }
            reads = self.run_process("upload-fastq-paired", inputs)
            ref_seq = self.prepare_ref_seq(
                fn="hs_b37_chr2_small.fasta.gz", species="Homo sapiens", build="b37",
            )
            bwa_index = self.run_process("bwa-index", {"ref_seq": ref_seq.id})
            master_file = self.prepare_amplicon_master_file()

        inputs = {
            "master_file": master_file.id,
            "genome": bwa_index.id,
            "reads": reads.id,
        }
        bwa_trim = self.run_process("align-bwa-trim", inputs)

        self.assertFile(bwa_trim, "stats", "bwa_trim_stats.txt")
        self.assertFile(bwa_trim, "bigwig", "bwa_trim_bigwig.bw")
        self.assertFields(bwa_trim, "species", "Homo sapiens")
        self.assertFields(bwa_trim, "build", "b37")

    @skipUnlessLargeFiles("56GSID_10k_mate1_RG.bam")
    @tag_process("amplicon-table")
    def test_amplicon_table(self):
        with self.preparation_stage():
            bam = self.run_process(
                "upload-bam",
                {
                    "src": join("large", "56GSID_10k_mate1_RG.bam"),
                    "species": "Homo sapiens",
                    "build": "b37",
                },
            )
            master_file = self.prepare_amplicon_master_file()

            coverage = self.run_process(
                "coveragebed", {"alignment": bam.id, "master_file": master_file.id,}
            )

            inputs = {
                "annotation": "56GSID.lf.finalvars.txt",
                "summary": "56GSID_1k.gatkHC_snpEff_summary.html",
                "snpeff_genes": "56GSID_1k.gatkHC_snpEff_genes.txt",
            }

            annot_variants = self.run_process("upload-snpeff", inputs)

        amplicon_table_inputs = {
            "master_file": master_file.id,
            "coverage": coverage.id,
            "annot_vars": [annot_variants.id],
        }

        table = self.run_process("amplicon-table", amplicon_table_inputs)
        self.assertJSON(
            table, table.output["variant_table"], "", "amplicon_table_output.json.gz"
        )

        amplicon_table_inputs = {
            "master_file": master_file.id,
            "coverage": coverage.id,
            "annot_vars": [annot_variants.id],
            "all_amplicons": True,
            "table_name": "All amplicons",
        }

        table = self.run_process("amplicon-table", amplicon_table_inputs)
        self.assertJSON(
            table,
            table.output["variant_table"],
            "",
            "amplicon_table_output_all.json.gz",
        )
