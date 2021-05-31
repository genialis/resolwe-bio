from pathlib import Path

from resolwe.flow.models import Process
from resolwe.test import tag_process

from resolwe_bio.utils.filter import filter_vcf_variable
from resolwe_bio.utils.test import BioProcessTestCase


class WgsProcessorTestCase(BioProcessTestCase):
    @tag_process("wgs-preprocess")
    def test_wgs_preprocess(self):
        def filter_startedon(line):
            """Filter stared on header line."""
            return line.startswith(b"# Started on:") or line.startswith(
                b"# MarkDuplicates"
            )

        base = Path("wgs")
        inputs = base / "input"
        outputs = base / "output"
        with self.preparation_stage():
            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": inputs / "hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )
            bwa_index = self.run_process("bwa-index", {"ref_seq": ref_seq.id})

            reads = self.prepare_paired_reads(
                mate1=[inputs / "TP53_1.fastq.gz"],
                mate2=[inputs / "TP53_2.fastq.gz"],
            )

            dbsnp = self.run_process(
                "upload-variants-vcf",
                {
                    "src": inputs / "dbsnp_TP53.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

        analysis_ready_bam = self.run_process(
            "wgs-preprocess",
            {
                "reads": reads.id,
                "ref_seq": ref_seq.id,
                "bwa_index": bwa_index.id,
                "known_sites": [dbsnp.id],
            },
        )

        self.assertFile(
            analysis_ready_bam,
            "stats",
            outputs / "wgs_preprocess_bam_stats.txt",
        )

        self.assertFile(
            analysis_ready_bam,
            "metrics_file",
            outputs / "wgs_preprocess_markdups_metrics.txt",
            file_filter=filter_startedon,
        )

    @tag_process("gatk-haplotypecaller-gvcf")
    def test_gatk_hc_gvcf(self):
        base = Path("wgs")
        inputs = base / "input"
        outputs = base / "output"
        with self.preparation_stage():
            input_bam = self.run_process(
                "upload-bam",
                {
                    "src": inputs / "analysis_ready.bam",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": inputs / "hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            intervals = self.run_process(
                "upload-bed",
                {
                    "src": inputs / "hg38.intervals.bed",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

        variants = self.run_process(
            "gatk-haplotypecaller-gvcf",
            {
                "bam": input_bam.id,
                "ref_seq": ref_seq.id,
                "options": {"intervals": intervals.id},
            },
        )

        self.assertFile(
            variants,
            "vcf",
            outputs / "variants.g.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )
        self.assertFields(variants, "build", "custom_build")
        self.assertFields(variants, "species", "Homo sapiens")

    @tag_process("gatk-genotype-gvcfs")
    def test_gatk_genotypegvcfs(self):
        base = Path("wgs")
        inputs = base / "input"
        outputs = base / "output"
        with self.preparation_stage():
            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": inputs / "hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            intervals = self.run_process(
                "upload-bed",
                {
                    "src": inputs / "hg38.intervals.bed",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            dbsnp = self.run_process(
                "upload-variants-vcf",
                {
                    "src": inputs / "dbsnp_TP53.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            # Mock upload gvcf process
            process = Process.objects.create(
                name="Upload GVCF mock process",
                requirements={
                    "expression-engine": "jinja",
                    "resources": {
                        "network": True,
                    },
                    "executor": {
                        "docker": {
                            "image": "public.ecr.aws/s4q6j6e8/resolwebio/base:ubuntu-20.04-03042021",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:variants:gvcf:",
                entity_type="sample",
                entity_descriptor_schema="sample",
                data_name="{{ gvcf.file }}",
                input_schema=[
                    {
                        "name": "gvcf",
                        "type": "basic:file:",
                    },
                    {
                        "name": "tabix",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "vcf",
                        "type": "basic:file:",
                    },
                    {
                        "name": "tbi",
                        "type": "basic:file:",
                    },
                    {
                        "name": "species",
                        "type": "basic:string:",
                    },
                    {
                        "name": "build",
                        "type": "basic:string:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ gvcf.file_temp|default(gvcf.file) }} {{ gvcf.file }} "g.vcf" "g.vcf" 0.1 compress
re-save-file vcf "${NAME}.g.vcf.gz"
re-import {{ tabix.file_temp|default(tabix.file) }} {{ tabix.file }} "g.vcf.gz.tbi" "g.vcf.gz.tbi" 0.1 extract
re-save-file tbi "${NAME}.g.vcf.gz.tbi"
re-save species "Homo sapiens"
re-save build "custom_build"
""",
                },
            )

            gvcf_input = {
                "gvcf": inputs / "variants.g.vcf.gz",
                "tabix": inputs / "variants.g.vcf.gz.tbi",
            }
            gvcf_1 = self.run_process(process.slug, gvcf_input)

            gvcf_input = {
                "gvcf": inputs / "variants2.g.vcf.gz",
                "tabix": inputs / "variants2.g.vcf.gz.tbi",
            }
            gvcf_2 = self.run_process(process.slug, gvcf_input)

        joint_variants = self.run_process(
            "gatk-genotype-gvcfs",
            {
                "gvcfs": [gvcf_1.id, gvcf_2.id],
                "ref_seq": ref_seq.id,
                "intervals": intervals.id,
                "dbsnp": dbsnp.id,
            },
        )

        self.assertFile(
            joint_variants,
            "vcf",
            outputs / "cohort_variants.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )
        self.assertFields(joint_variants, "build", "custom_build")
        self.assertFields(joint_variants, "species", "Homo sapiens")
