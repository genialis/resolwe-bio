from pathlib import Path

from resolwe.flow.models import Process
from resolwe.test import tag_process

from resolwe_bio.utils.filter import filter_vcf_variable
from resolwe_bio.utils.test import BioProcessTestCase


class WgsProcessorTestCase(BioProcessTestCase):
    @tag_process("wgs-preprocess-bwa2")
    def test_wgs_preprocess_bwa2(self):
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
            bwa_index = self.run_process("bwamem2-index", {"ref_seq": ref_seq.id})

            reads = self.prepare_paired_reads(
                mate1=[inputs / "TP53_1.fastq.gz"],
                mate2=[inputs / "TP53_2.fastq.gz"],
            )

            aligned_reads = self.run_process(
                "upload-bam",
                {
                    "src": inputs / "TP53_aligned_reads.bam",
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

        analysis_ready_bam = self.run_process(
            "wgs-preprocess-bwa2",
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

        # test the process with aligned BAM provided as an input
        analysis_ready_bam = self.run_process(
            "wgs-preprocess-bwa2",
            {
                "aligned_reads": aligned_reads.id,
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

    @tag_process("gatk-genomicsdb-import")
    def test_gatk_genomicsdb(self):
        base = Path("wgs")
        inputs = base / "input"
        with self.preparation_stage():
            intervals = self.run_process(
                "upload-bed",
                {
                    "src": inputs / "hg38.intervals.bed",
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
                            "image": "public.ecr.aws/genialis/resolwebio/base:ubuntu-22.04-14112023",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:variants:gvcf:",
                entity_type="sample",
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

        database = self.run_process(
            "gatk-genomicsdb-import",
            {
                "gvcfs": [gvcf_1.id],
                "intervals": intervals.id,
            },
        )

        self.assertFields(database, "build", "custom_build")
        self.assertFields(database, "species", "Homo sapiens")

        database2 = self.run_process(
            "gatk-genomicsdb-import",
            {
                "gvcfs": [gvcf_2.id],
                "use_existing": True,
                "existing_db": database.id,
            },
        )
        self.assertFields(database2, "build", "custom_build")
        self.assertFields(database2, "species", "Homo sapiens")

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
                            "image": "public.ecr.aws/genialis/resolwebio/base:ubuntu-22.04-14112023",
                        },
                    },
                },
                contributor=self.contributor,
                type="data:variants:gvcf:",
                entity_type="sample",
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

            database = self.run_process(
                "gatk-genomicsdb-import",
                {
                    "gvcfs": [gvcf_1.id, gvcf_2.id],
                    "intervals": intervals.id,
                },
            )

        joint_variants = self.run_process(
            "gatk-genotype-gvcfs",
            {
                "database": database.id,
                "ref_seq": ref_seq.id,
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

    @tag_process("gatk-merge-vcfs")
    def test_gatk_merge_vcfs(self):
        base = Path("wgs")
        inputs = base / "input"
        outputs = base / "output"
        with self.preparation_stage():
            vcf_1 = self.run_process(
                "upload-variants-vcf",
                {
                    "src": inputs / "vcf_1.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            vcf_2 = self.run_process(
                "upload-variants-vcf",
                {
                    "src": inputs / "vcf_2.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

        merged_vcfs = self.run_process(
            "gatk-merge-vcfs",
            {"vcfs": [vcf_1.id, vcf_2.id]},
        )

        self.assertFile(
            merged_vcfs,
            "vcf",
            outputs / "combined_variants.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )
        self.assertFields(merged_vcfs, "build", "custom_build")
        self.assertFields(merged_vcfs, "species", "Homo sapiens")

    @tag_process("gatk-select-variants", "gatk-select-variants-single")
    def test_gatk_select_variants(self):
        base = Path("wgs")
        inputs = base / "input"
        outputs = base / "output"
        with self.preparation_stage():
            vcf = self.run_process(
                "upload-variants-vcf",
                {
                    "src": inputs / "vcf_1.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            intervals = self.run_process(
                "upload-bed",
                {
                    "src": inputs / "chr17_single_interval.bed",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

            vcf_filtered = self.run_process(
                "upload-variants-vcf",
                {
                    "src": inputs / "filtered_variants.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

        selected_variants = self.run_process(
            process_slug="gatk-select-variants",
            input_={
                "vcf": vcf.id,
                "intervals": intervals.id,
                "select_type": ["SNP", "INDEL"],
            },
        )

        self.assertFile(
            selected_variants,
            "vcf",
            outputs / "selected_variants.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )
        self.assertFields(selected_variants, "build", "custom_build")
        self.assertFields(selected_variants, "species", "Homo sapiens")

        exclude_filtered = self.run_process(
            process_slug="gatk-select-variants-single",
            input_={
                "vcf": vcf_filtered.id,
                "select_type": ["SNP"],
                "exclude_filtered": True,
            },
        )

        self.assertFile(
            exclude_filtered,
            "vcf",
            outputs / "selected_variants_filtered.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )
