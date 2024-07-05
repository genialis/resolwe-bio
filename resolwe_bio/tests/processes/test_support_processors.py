import os
from pathlib import Path

from resolwe.flow.models import Data, Process
from resolwe.flow.models.entity import Relation, RelationPartition, RelationType
from resolwe.permissions.models import Permission
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.models import Sample
from resolwe_bio.utils.filter import filter_comment_lines, filter_rnaseqc_metrics
from resolwe_bio.utils.test import KBBioProcessTestCase
from resolwe_bio.variants.models import Variant, VariantCall, VariantExperiment


class SupportProcessorTestCase(KBBioProcessTestCase):
    fixtures = ["relationtypes.yaml"]

    @tag_process("bam-split")
    def test_bam_split(self):
        bam_input = Path("test_bam", "input")
        bam_output = Path("test_bam", "output")

        with self.preparation_stage():
            bam = self.prepare_bam(
                fn=bam_input / "hybrid.bam", species="Mus musculus", build="mm10_dm6"
            )

            header = self.run_process(
                "upload-header-sam", {"src": bam_input / "mm10_header.sam"}
            )
            header2 = self.run_process(
                "upload-header-sam", {"src": bam_input / "dm6_header.sam"}
            )

        inputs = {
            "bam": bam.id,
        }
        bam1 = self.run_process("bam-split", inputs)
        bam2 = Data.objects.last()

        self.assertFile(bam1, "bam", bam_output / "hybrid_mm10.bam")
        self.assertFile(bam1, "bai", bam_output / "hybrid_mm10.bam.bai")
        self.assertFields(bam1, "species", "Mus musculus")
        self.assertFields(bam1, "build", "mm10")
        self.assertFile(bam2, "bam", bam_output / "hybrid_dm6.bam")
        self.assertFile(bam2, "bai", bam_output / "hybrid_dm6.bam.bai")
        self.assertFields(bam2, "species", "Drosophila melanogaster")
        self.assertFields(bam2, "build", "dm6")

        inputs["header"] = header.id
        inputs["header2"] = header2.id
        bam1 = self.run_process("bam-split", inputs)
        bam2 = Data.objects.last()

        self.assertFile(bam1, "bam", bam_output / "hybrid_mm10.bam")
        self.assertFile(bam1, "bai", bam_output / "hybrid_mm10.bam.bai")
        self.assertFields(bam1, "species", "Mus musculus")
        self.assertFields(bam1, "build", "mm10")
        self.assertFile(bam2, "bam", bam_output / "hybrid_dm6.bam")
        self.assertFile(bam2, "bai", bam_output / "hybrid_dm6.bam.bai")
        self.assertFields(bam2, "species", "Drosophila melanogaster")
        self.assertFields(bam2, "build", "dm6")

    @tag_process("gff-to-gtf")
    def test_gff_to_gtf(self):
        with self.preparation_stage():
            annotation = self.prepare_annotation_gff()

        gff_to_gtf = self.run_process("gff-to-gtf", {"annotation": annotation.id})
        self.assertFile(gff_to_gtf, "annot", "gff_to_gtf_annotation.gtf")
        del gff_to_gtf.output["annot_sorted_track_jbrowse"][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFields(
            gff_to_gtf,
            "annot_sorted_track_jbrowse",
            {"refs": ["tracks/annotation"], "file": "trackList.json"},
        )

    @with_resolwe_host
    @tag_process("archive-samples")
    def test_ars(self):
        with self.preparation_stage():
            txt_file = self.run_process("upload-file", {"src": "markdup_stats.txt"})
            bam_input = {
                "src": "bamplot_alignment.bam",
                "species": "Mus musculus",
                "build": "GRCh38 _ens90",
            }
            bam = self.run_process("upload-bam", bam_input)

            read_inputs = {"src": ["rRNA forw.fastq.gz", "rRNA_rew.fastq.gz"]}
            reads = self.run_process("upload-fastq-single", read_inputs)

            vcf_input = {
                "src": "igv_human.lf.vcf",
                "species": "Homo sapiens",
                "build": "b37",
            }
            vcf = self.run_process("upload-variants-vcf", vcf_input)

            expression_1 = self.prepare_expression(f_exp="exp_1_rc.tab.gz", f_type="RC")
            expression_2 = self.prepare_expression(f_exp="exp_2_rc.tab.gz", f_type="RC")
            expression_5 = self.prepare_expression(f_exp="exp_5_rc.tab.gz", f_type="RC")
            expression_3 = self.prepare_expression(
                f_rc="exp_2_rc.tab.gz", f_exp="exp_2_tpm.tab.gz", f_type="TPM"
            )
            expression_4 = self.prepare_expression(
                f_rc="exp_2_rc.tab.gz", f_exp="exp_2_tpm.tab.gz", f_type="TPM"
            )

            multiqc = self.run_process(
                "multiqc",
                {
                    "data": [
                        bam.id,
                        reads.id,
                    ],
                    "advanced": {
                        "dirs": True,
                        "config": True,
                    },
                },
            )

        self.run_process(
            "archive-samples",
            {
                "data": [
                    txt_file.id,
                    bam.id,
                    reads.id,
                    vcf.id,
                    expression_1.id,
                    expression_2.id,
                    expression_5.id,
                    expression_4.id,
                    expression_3.id,
                    multiqc.id,
                ],
                "fields": [
                    "file",
                    "bam",
                    "bai",
                    "fastq",
                    "fastqc_url",
                    "fastqc_archive",
                    "vcf",
                    "exp_set",
                    "report",
                ],
            },
        )

    @with_resolwe_host
    @tag_process("archive-samples")
    def test_archive_samples_exp_set(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(f_exp="exp_1_rc.tab.gz", f_type="RC")
            expression_2 = self.prepare_expression(
                f_exp="exp_2_tpm.tab.gz", f_type="TPM", f_rc="exp_2_rc.tab.gz"
            )

        inputs = {
            "data": [
                expression_1.pk,
                expression_2.pk,
            ],
            "fields": [
                "exp_set",
            ],
        }
        self.run_process("archive-samples", inputs)
        # Structured zip files are not supported by assertFile. When implemented, add here
        # self.assertFile(_, 'archive', 'test_archive_samples_exp_set.zip', compression='zip').

    @with_resolwe_host
    @tag_process("archive-samples")
    def test_archive_samples_exp(self):
        with self.preparation_stage():
            # Upload expression without exp_set output.
            Process.objects.create(
                slug="exp",
                type="data:expression:exp:",
                contributor=self.contributor,
                requirements={"expression-engine": "jinja"},
                data_name="Upload expression into exp output field",
                entity_type="sample",
                input_schema=[
                    {
                        "name": "exp",
                        "type": "basic:file:",
                    }
                ],
                output_schema=[
                    {
                        "name": "exp",
                        "type": "basic:file:",
                    },
                    {
                        "name": "exp_type",
                        "type": "basic:string:",
                    },
                    {
                        "name": "build",
                        "type": "basic:string:",
                    },
                    {
                        "name": "species",
                        "type": "basic:string:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": """
                        re-import {{ exp.file_temp|default(exp.file) }} {{ exp.file }} "tab|gz" tab 1.0 compress
                        re-save-file exp "${NAME}.tab.gz"
                        re-save exp_type RC
                        re-save build ens
                        re-save species "Homo sapiens"
                    """,
                },
            )
            inputs = {
                "exp": "exp_1_rc.tab.gz",
            }
            expression_1 = self.run_process("exp", inputs)
            expression_2 = self.prepare_expression(
                f_exp="exp_2_tpm.tab.gz", f_type="TPM", f_rc="exp_2_rc.tab.gz"
            )

        inputs = {
            "data": [
                expression_1.pk,
                expression_2.pk,
            ],
            "fields": [
                "exp_set",
            ],
        }
        self.run_process("archive-samples", inputs)
        # Structured zip files are not supported by assertFile. When implemented, add here
        # self.assertFile(_, 'archive', 'test_archive_samples_exp.zip', compression='zip').

    @tag_process("prepare-geo-chipseq")
    def test_prepare_geo_chipseq(self):
        with self.preparation_stage():
            reads_1 = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )
            reads_2 = self.prepare_reads()
            reads_3 = self.prepare_reads(["SRR2124780_1 1k.fastq.gz"])

            macs14_case_bam = self.prepare_bam(
                fn="macs14_case.bam", species="Homo sapiens", build="hg19"
            )
            macs14_control_bam = self.prepare_bam(
                fn="macs14_control.bam", species="Homo sapiens", build="hg19"
            )

            macs2_case_bam = self.prepare_bam(
                fn="macs2/input/SRR5675973_chr17.bam",
                species="Homo sapiens",
                build="hg19",
            )
            macs2_control_bam = self.prepare_bam(
                fn="macs2/input/SRR5675974_chr17.bam",
                species="Homo sapiens",
                build="hg19",
            )

            # Run macs14
            inputs = {"treatment": macs14_case_bam.id, "control": macs14_control_bam.id}
            macs14_1 = self.run_process("macs14", inputs)
            macs14_1.entity = reads_1.entity
            macs14_1.save()

            macs14_control_bam.entity = reads_2.entity
            macs14_control_bam.save()

            # Run macs14 without control/background sample
            del inputs["control"]
            macs14_2 = self.run_process("macs14", inputs)
            macs14_2.entity = reads_3.entity
            macs14_2.save()

            # Run macs2
            inputs = {
                "case": macs2_case_bam.id,
                "control": macs2_control_bam.id,
                "settings": {
                    "extsize": 298,
                    "nomodel": True,
                    "bedgraph": True,
                },
            }
            macs2_1 = self.run_process("macs2-callpeak", inputs)
            macs2_1.entity = reads_1.entity
            macs2_1.save()

            macs2_control_bam.entity = reads_2.entity
            macs2_control_bam.save()

            # Run macs2 without control/background sample
            del inputs["control"]
            macs2_2 = self.run_process("macs2-callpeak", inputs)
            macs2_2.entity = reads_3.entity
            macs2_2.save()

        inputs = {
            "reads": [reads_1.id, reads_2.id, reads_3.id],
            "macs": [macs14_1.id, macs14_2.id],
            "name": "prepare_geo",
        }
        prepare_geo_chipseq = self.run_process("prepare-geo-chipseq", inputs)

        self.assertFile(prepare_geo_chipseq, "table", "prepare_geo_ChIP-Seq_macs14.txt")

        inputs = {
            "reads": [reads_1.id, reads_2.id, reads_3.id],
            "macs": [macs2_1.id, macs2_2.id],
            "name": "prepare_geo",
        }
        prepare_geo_chipseq = self.run_process("prepare-geo-chipseq", inputs)

        self.assertFile(prepare_geo_chipseq, "table", "prepare_geo_ChIP-Seq_macs2.txt")

    @with_resolwe_host
    @tag_process("prepare-geo-rnaseq")
    def test_prepare_geo_rnaseq(self):
        with self.preparation_stage():
            reads_1 = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )
            reads_2 = self.prepare_reads()

            expression_1 = self.prepare_expression(
                f_rc="exp_1_rc.tab.gz", f_exp="exp_1_tpm.tab.gz", f_type="TPM"
            )
            expression_2 = self.prepare_expression(
                f_rc="exp_2_rc.tab.gz", f_exp="exp_2_tpm.tab.gz", f_type="TPM"
            )

            # Add expressions to reads samples
            expression_1.entity = reads_1.entity
            expression_1.save()
            expression_2.entity = reads_2.entity
            expression_2.save()

        inputs = {
            "reads": [reads_1.id, reads_2.id],
            "expressions": [expression_1.pk, expression_2.pk],
            "name": "prepare_geo",
        }
        prepare_geo_rnaseq = self.run_process("prepare-geo-rnaseq", inputs)

        self.assertFile(prepare_geo_rnaseq, "table", "prepare_geo_RNA-Seq.txt")

    @tag_process("library-strandedness")
    def test_library_strandedness(self):
        with self.preparation_stage():
            cds = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "salmon_cds.fa.gz",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            inputs = {
                "nucl": cds.id,
                "source": "ENSEMBL",
                "species": "Homo sapiens",
                "build": "ens_90",
            }
            salmon_index = self.run_process("salmon-index", inputs)

            single_reads = self.prepare_reads(["reads rsem.fq.gz"])
            paired_reads = self.prepare_paired_reads(
                mate1=["reads rsem.fq.gz"], mate2=["reads rsem2.fq.gz"]
            )

        single_input = {
            "reads": single_reads.id,
            "salmon_index": salmon_index.id,
        }

        lib_strandedness_single = self.run_process("library-strandedness", single_input)
        self.assertFields(lib_strandedness_single, "strandedness", "U")
        self.assertFields(lib_strandedness_single, "fragment_ratio", 1.0)

        paired_input = {
            "reads": paired_reads.id,
            "salmon_index": salmon_index.id,
        }

        lib_strandedness_paired = self.run_process("library-strandedness", paired_input)
        self.assertFields(lib_strandedness_paired, "strandedness", "IU")
        self.assertFields(lib_strandedness_paired, "fragment_ratio", 1.0)

    @tag_process("multiqc")
    def test_multiqc_duplicate_inputs(self):
        with self.preparation_stage():
            reads = self.run_processor(
                "upload-fastq-single",
                {"src": ["chr1_single.fastq.gz"]},
            )

            annotation = self.run_process(
                "upload-gtf",
                {
                    "src": "chr1_region.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            genome_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "chr1_region.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )
            star_index = self.run_process(
                "alignment-star-index",
                {
                    "annotation": annotation.id,
                    "ref_seq": genome_fasta.id,
                },
            )

            star_alignment_1 = self.run_process(
                "alignment-star",
                {
                    "genome": star_index.id,
                    "reads": reads.id,
                    "gene_counts": True,
                },
            )

            star_alignment_2 = self.run_process(
                "alignment-star",
                {
                    "genome": star_index.id,
                    "reads": reads.id,
                    "gene_counts": True,
                },
            )

            multiqc = self.run_process(
                "multiqc",
                {
                    "data": [
                        reads.id,
                        star_alignment_1.id,
                        star_alignment_2.id,
                    ],
                    "advanced": {
                        "dirs": True,
                        "config": True,
                        "cl_config": "fn_ignore_dirs: [tmp]",
                    },
                },
            )
        self.assertFileExists(multiqc, "report")

    @tag_process("multiqc")
    def test_multiqc(self):
        with self.preparation_stage():
            reads = self.run_processor(
                "upload-fastq-single",
                {"src": ["chr1_single.fastq.gz"]},
            )

            paired_reads = self.prepare_paired_reads(
                ["chr1_paired_R1.fastq.gz"],
                ["chr1_paired_R2.fastq.gz"],
            )

            filtered_reads = self.run_process(
                "bbduk-paired", {"reads": paired_reads.id}
            )

            bam_samtools = self.run_process(
                "upload-bam-indexed",
                {
                    "src": Path("test_bam", "output", "alignment_position_sorted.bam"),
                    "src2": Path(
                        "test_bam", "output", "alignment_position_sorted.bam.bai"
                    ),
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

            annotation = self.run_process(
                "upload-gtf",
                {
                    "src": "chr1_region.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            genome_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "chr1_region.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )
            star_index = self.run_process(
                "alignment-star-index",
                {
                    "annotation": annotation.id,
                    "ref_seq": genome_fasta.id,
                },
            )

            star_alignment = self.run_process(
                "alignment-star",
                {
                    "genome": star_index.id,
                    "reads": paired_reads.id,
                    "gene_counts": True,
                },
            )

            star_quantification = self.run_process(
                "star-quantification",
                {
                    "aligned_reads": star_alignment.id,
                    "annotation": annotation.id,
                },
            )

            samtools_idxstats = self.run_process(
                "samtools-idxstats",
                {
                    "alignment": star_alignment.id,
                },
            )

            qorts_report = self.run_process(
                "qorts-qc",
                {
                    "alignment": star_alignment.id,
                    "annotation": annotation.id,
                    "options": {
                        "maxPhredScore": 42,
                    },
                },
            )

            # BED file is not part of a sample entity. Test if MultiQC process
            # correctly skips this input data object
            bed = self.run_process(
                "upload-bed",
                {"src": "good.bed", "species": "Homo sapiens", "build": "hg19"},
            )

        multiqc = self.run_process(
            "multiqc",
            {
                "data": [
                    reads.id,
                    paired_reads.id,
                    filtered_reads.id,
                    bam_samtools.id,
                    star_alignment.id,
                    star_quantification.id,
                    samtools_idxstats.id,
                    qorts_report.id,
                    bed.id,
                ],
                "advanced": {
                    "dirs": True,
                    "config": True,
                    "cl_config": "fn_ignore_dirs: [tmp]",
                },
            },
        )
        self.assertFileExists(multiqc, "report")

    @tag_process("multiqc")
    def test_multiqc_wgbs(self):
        with self.preparation_stage():

            def set_sample_name(data, sample_name):
                """Set sample name."""
                sample = Sample.objects.get(data=data)
                sample.name = sample_name
                sample.save()

            process = Process.objects.create(
                name="Upload bsrate data mock process",
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
                entity_type="sample",
                contributor=self.contributor,
                type="data:wgbs:bsrate:",
                input_schema=[
                    {
                        "name": "src",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "report",
                        "type": "basic:file:",
                    }
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "txt" "txt" 0.1 extract
re-save-file report "${NAME}".txt
""",
                },
            )

            bsrate = self.run_process(
                process.slug,
                {
                    "src": os.path.join(
                        "wgbs", "output", "Escherichia_phage_bsrate.txt"
                    ),
                },
            )
            set_sample_name(bsrate, "Bsrate test")

        multiqc = self.run_process(
            "multiqc",
            {
                "data": [
                    bsrate.id,
                ],
                "advanced": {
                    "dirs": True,
                    "config": True,
                },
            },
        )
        self.assertFileExists(multiqc, "report")

    @tag_process("multiqc")
    def test_multiqc_markdup(self):
        with self.preparation_stage():

            def set_sample_name(data, sample_name):
                """Set sample name."""
                sample = Sample.objects.get(data=data)
                sample.name = sample_name
                sample.save()

            process = Process.objects.create(
                name="Upload walt data mock process",
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
                entity_type="sample",
                contributor=self.contributor,
                type="data:alignment:bam:walt:",
                input_schema=[
                    {
                        "name": "src",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "duplicates_report",
                        "type": "basic:file:",
                    }
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "txt" "txt" 0.1 extract
re-save-file duplicates_report "${NAME}".txt
""",
                },
            )

            walt = self.run_process(
                process.slug,
                {"src": "markdup_stats.txt"},
            )
            set_sample_name(walt, "Walt markdup test")

        multiqc = self.run_process(
            "multiqc",
            {
                "data": [
                    walt.id,
                ],
                "advanced": {
                    "dirs": True,
                    "config": True,
                },
            },
        )
        self.assertFileExists(multiqc, "report")

    @tag_process("multiqc")
    def test_multiqc_chipqc(self):
        with self.preparation_stage():

            def set_sample_name(data, sample_name):
                """Set sample name."""
                sample = Sample.objects.get(data=data)
                sample.name = sample_name
                sample.save()

            process = Process.objects.create(
                name="Upload chipqc data mock process",
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
                entity_type="sample",
                contributor=self.contributor,
                type="data:chipqc:",
                input_schema=[
                    {
                        "name": "src",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "ccplot",
                        "type": "basic:file:",
                    },
                    {
                        "name": "coverage_histogram",
                        "type": "basic:file:",
                    },
                    {
                        "name": "peak_profile",
                        "type": "basic:file:",
                    },
                    {
                        "name": "peaks_barplot",
                        "type": "basic:file:",
                    },
                    {
                        "name": "peaks_density_plot",
                        "type": "basic:file:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "png" "png" 0.1 extract
cp "${NAME}".png CCPlot_mqc.png
cp "${NAME}".png CoverageHistogramPlot_mqc.png
cp "${NAME}".png Rip_mqc.png
cp "${NAME}".png Rap_mqc.png
re-save-file ccplot CCPlot_mqc.png
re-save-file coverage_histogram CoverageHistogramPlot_mqc.png
re-save-file peak_profile "${NAME}".png
re-save-file peaks_barplot Rip_mqc.png
re-save-file peaks_density_plot Rap_mqc.png
""",
                },
            )

            peak_qc = Process.objects.create(
                name="Upload Post-Peak QC mock process",
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
                entity_type="sample",
                contributor=self.contributor,
                type="data:chipseq:callpeak:macs2:",
                input_schema=[
                    {
                        "name": "src",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "chip_qc",
                        "type": "basic:file:",
                    },
                    {
                        "name": "called_peaks",
                        "type": "basic:file:",
                    },
                    {
                        "name": "case_prepeak_qc",
                        "type": "basic:file:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "txt" "txt" 0.1 extract
re-save-file chip_qc "${NAME}".txt
re-save-file called_peaks "${NAME}".txt
re-save-file case_prepeak_qc "${NAME}".txt
""",
                },
            )

            chipqc = self.run_process(
                process.slug,
                {
                    "src": os.path.join("chipqc", "output", "PeakProfile_mqc.png"),
                },
            )
            set_sample_name(chipqc, "ChipQC test.fastq.gz")

            postpeak_qc_report = self.run_process(
                peak_qc.slug,
                {"src": os.path.join("chipqc", "input", "postpeak_qc_report.txt")},
            )
            set_sample_name(postpeak_qc_report, "ChipQC test.fastq.gz")

        multiqc = self.run_process(
            "multiqc",
            {
                "data": [
                    chipqc.id,
                    postpeak_qc_report.id,
                ],
                "advanced": {
                    "dirs": True,
                    "config": True,
                },
            },
        )
        self.assertFileExists(multiqc, "report")

    @tag_process("multiqc")
    def test_multiqc_nanostring(self):
        with self.preparation_stage():

            def set_sample_name(data, sample_name):
                """Set sample name."""
                sample = Sample.objects.get(data=data)
                sample.name = sample_name
                sample.save()

            process = Process.objects.create(
                name="Upload rcc data mock process",
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
                entity_type="sample",
                contributor=self.contributor,
                type="data:nanostring:rcc:",
                input_schema=[
                    {
                        "name": "src",
                        "type": "basic:file:",
                    },
                    {
                        "name": "attr",
                        "type": "basic:file:",
                    },
                ],
                output_schema=[
                    {
                        "name": "sample_qc",
                        "type": "basic:file:",
                    },
                    {
                        "name": "lane_attributes",
                        "type": "basic:file:",
                    },
                ],
                run={
                    "language": "bash",
                    "program": r"""
re-import {{ src.file_temp|default(src.file) }} {{ src.file }} "txt" "txt" 0.1 extract
re-save-file sample_qc "${NAME}".txt
re-import {{ attr.file_temp|default(attr.file) }} {{ attr.file }} "txt" "txt" 0.1 extract
re-save-file lane_attributes "${NAME}".txt
""",
                },
            )

            rcc_1 = self.run_process(
                process.slug,
                {
                    "src": "nanostring_qc.txt",
                    "attr": "nanostring_lane_attr.txt",
                },
            )
            set_sample_name(rcc_1, "sample_1")

            rcc_2 = self.run_process(
                process.slug,
                {
                    "src": "nanostring_qc_2.txt",
                    "attr": "nanostring_lane_attr.txt",
                },
            )
            set_sample_name(rcc_2, "sample_with_NA")

        multiqc = self.run_process(
            "multiqc",
            {
                "data": [rcc_1.id, rcc_2.id],
                "advanced": {
                    "dirs": True,
                    "config": True,
                },
            },
        )
        self.assertFileExists(multiqc, "report")

    @tag_process("seqtk-sample-single", "seqtk-sample-paired")
    def test_seqtk_sample(self):
        input_folder = Path("seqtk") / "input"
        output_folder = Path("seqtk") / "output"
        with self.preparation_stage():
            reads = self.run_processor(
                "upload-fastq-single",
                {
                    "src": [
                        input_folder
                        / "hs_single_bbduk_star_htseq_reads_single.fastq.gz"
                    ]
                },
            )

            paired_reads = self.prepare_paired_reads(
                [input_folder / "hs_paired_R1_workflow_bbduk_star_htseq.fastq.gz"],
                [input_folder / "hs_paired_R2_workflow_bbduk_star_htseq.fastq.gz"],
            )

        inputs_single = {
            "reads": reads.id,
            "n_reads": 42,
            "advanced": {
                "seed": 42,
            },
        }

        seqtk_single = self.run_process("seqtk-sample-single", inputs_single)
        self.assertFiles(
            seqtk_single,
            "fastq",
            [output_folder / "seqtk_subsampled_reads_single_end.fastq.gz"],
            compression="gzip",
        )

        inputs_paired = {
            "reads": paired_reads.id,
            "advanced": {
                "seed": 42,
                "fraction": 0.25,
            },
        }

        seqtk_paired = self.run_process("seqtk-sample-paired", inputs_paired)
        self.assertFiles(
            seqtk_paired,
            "fastq",
            [output_folder / "seqtk_subsampled_reads_paired_end_mate1.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            seqtk_paired,
            "fastq2",
            [output_folder / "seqtk_subsampled_reads_paired_end_mate2.fastq.gz"],
            compression="gzip",
        )

    @with_resolwe_host
    @tag_process("spikein-qc")
    def test_spikein_pairwise(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(
                f_exp="exp1_cpm_ercc_sirv.tab.gz",
                f_type="CPM",
                name="Sample 1",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_2 = self.prepare_expression(
                f_exp="exp2_cpm_ercc_sirv.tab.gz",
                f_type="CPM",
                name="Sample 2",
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_3 = self.prepare_expression(
                f_exp="exp3_cpm_ercc_sirv.tab.gz",
                f_type="CPM",
                name="Sample3.txt",  # Test for sample names that might look like filename extensions
                source="ENSEMBL",
                species="Homo sapiens",
            )
            expression_4 = self.prepare_expression(
                f_exp="exp4_cpm_ercc_sirv.tab.gz",
                f_type="CPM",
                name="Sample without ERCC",
                source="ENSEMBL",
                species="Homo sapiens",
            )

        # SIRV Set 3
        sirv_set3 = self.run_process(
            "spikein-qc",
            {
                "samples": [
                    expression_1.pk,
                    expression_2.pk,
                    expression_3.pk,
                    expression_4.pk,
                ],
                "mix": "sirv_set3",
            },
        )

        self.assertEqual(
            sirv_set3.process_warning,
            [
                "All ERCC spike-ins have zero expression in sample Sample without ERCC",
            ],
        )

        self.assertFilesExist(sirv_set3, "plots")
        expected_names = [item["file"] for item in sirv_set3.output["plots"]]
        self.assertEqual(
            expected_names,
            [
                "Sample 1 (ERCC spike-in's).png",
                "Sample 2 (ERCC spike-in's).png",
                "Sample3.txt (ERCC spike-in's).png",
            ],
        )

        self.assertFileExists(sirv_set3, "report")

        self.assertFileExists(sirv_set3, "report_zip")

    @tag_process("rnaseqc-qc")
    def test_rnaseqc_qc(self):
        base = Path("rnaseqc")
        outputs = base / "output"
        with self.preparation_stage():
            alignment_ensembl = self.run_process(
                "upload-bam",
                {
                    "src": "chr1_region.bam",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )
            alignment_ucsc = self.run_process(
                "upload-bam",
                {
                    "src": "chr1_region_ucsc.bam",
                    "species": "Homo sapiens",
                    "build": "hg_19",
                },
            )

            cds = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "chr1_region.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )
            inputs_salmon_index = {
                "nucl": cds.id,
                "source": "ENSEMBL",
                "species": "Homo sapiens",
                "build": "ens_90",
            }
            salmon_index = self.run_process("salmon-index", inputs_salmon_index)

            annotation_ensembl = self.run_process(
                "upload-gtf",
                {
                    "src": "chr1_region.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            annotation_ucsc = self.run_process(
                "upload-gtf",
                {
                    "src": "chr1_region_UCSC.gtf.gz",
                    "source": "UCSC",
                    "species": "Homo sapiens",
                    "build": "hg_19",
                },
            )

        inputs_ensembl = {
            "annotation": annotation_ensembl.id,
            "alignment": alignment_ensembl.id,
            "strand_detection_options": {
                "stranded": "auto",
                "cdna_index": salmon_index.id,
            },
        }

        rnaseqc_report_ensembl = self.run_process("rnaseqc-qc", inputs_ensembl)

        self.assertFile(
            rnaseqc_report_ensembl,
            "metrics",
            outputs / "chr1_region.bam.metrics_ENSEMBL.tsv",
        )

        inputs_ucsc = {
            "annotation": annotation_ucsc.id,
            "alignment": alignment_ucsc.id,
            "rnaseqc_options": {"detection_threshold": 2},
        }

        rnaseqc_report_ucsc = self.run_process("rnaseqc-qc", inputs_ucsc)

        self.assertFile(
            rnaseqc_report_ucsc,
            "metrics",
            outputs / "chr1_region_ucsc.bam.metrics_UCSC.tsv",
            # Because testing on Mac and Linux yields nan and -nan, respectively,
            # we are, for now, omitting this parameter from tests.
            file_filter=filter_rnaseqc_metrics,
        )

    @tag_process("qorts-qc")
    def test_qorts_qc(self):
        with self.preparation_stage():
            alignment = self.run_process(
                "upload-bam",
                {
                    "src": "qorts/input/hs paired.bam",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            cds = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "qorts/input/salmon_cds.fa.gz",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )
            inputs = {
                "nucl": cds.id,
                "source": "ENSEMBL",
                "species": "Homo sapiens",
                "build": "ens_90",
            }
            salmon_index = self.run_process("salmon-index", inputs)

            annotation = self.run_process(
                "upload-gtf",
                {
                    "src": "qorts/input/hs annotation.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

        inputs = {
            "alignment": alignment.id,
            "annotation": annotation.id,
            "options": {
                "stranded": "auto",
                "cdna_index": salmon_index.id,
                "adjustPhredScore": 31,
            },
        }
        qorts_report = self.run_process("qorts-qc", inputs)
        self.assertFileExists(qorts_report, "plot")
        self.assertFileExists(qorts_report, "summary")
        self.assertFileExists(qorts_report, "qorts_data")

    @tag_process("samtools-idxstats")
    def test_samtools_idxstats(self):
        with self.preparation_stage():
            alignment = self.prepare_bam(
                fn=Path("test_bam", "output", "alignment_position_sorted.bam")
            )

        idxstats = self.run_process("samtools-idxstats", {"alignment": alignment.id})
        self.assertFile(idxstats, "report", "samtools_idxstats_report.txt")

    @tag_process("umi-tools-dedup")
    def test_umi_tools_dedup(self):
        with self.preparation_stage():
            bam_single = self.run_process(
                "upload-bam",
                {
                    "src": "./corall/input/corall_single.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            bam_paired = self.run_process(
                "upload-bam",
                {
                    "src": "./corall/input/corall_paired.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

        dedup_single = self.run_process("umi-tools-dedup", {"alignment": bam_single.id})
        self.assertFile(dedup_single, "stats", "./corall/output/dedup_single_stats.txt")

        dedup_paired = self.run_process("umi-tools-dedup", {"alignment": bam_paired.id})
        self.assertFile(dedup_paired, "stats", "./corall/output/dedup_paired_stats.txt")

    @tag_process("seqtk-rev-complement-single", "seqtk-rev-complement-paired")
    def test_reverse_complement(self):
        with self.preparation_stage():
            single_reads = self.prepare_reads(["hs_slamseq_R1.fastq.gz"])

            paired_reads = self.prepare_paired_reads(
                ["hs_slamseq_R1.fastq.gz"], ["hs_slamseq_R2.fastq.gz"]
            )

        revcomp_single = self.run_process(
            "seqtk-rev-complement-single", {"reads": single_reads.id}
        )

        self.assertFiles(
            revcomp_single,
            "fastq",
            ["hs_slamseq_R1_complemented.fastq.gz"],
            compression="gzip",
        )

        revcomp_paired = self.run_process(
            "seqtk-rev-complement-paired",
            {"reads": paired_reads.id, "select_mate": "Mate 1"},
        )

        self.assertFiles(
            revcomp_paired,
            "fastq",
            ["hs_slamseq_R1_complemented.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            revcomp_paired, "fastq2", ["hs_slamseq_R2.fastq.gz"], compression="gzip"
        )

    @tag_process("alignment-summary")
    def test_alignment_summary(self):
        with self.preparation_stage():
            bam = self.run_process(
                "upload-bam",
                {
                    "src": "bamclipper/output/TP53.primerclipped.bam",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )
            genome = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "bqsr/input/hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )
            adapters = self.prepare_ref_seq()

        alignment_summary = self.run_process(
            "alignment-summary",
            {
                "bam": bam.id,
                "genome": genome.id,
                "adapters": adapters.id,
            },
        )

        self.assertFile(
            alignment_summary,
            "report",
            "hs_gatk_alignment_summary_metrics.txt",
            file_filter=filter_comment_lines,
        )

    @tag_process("insert-size")
    def test_insert_size(self):
        with self.preparation_stage():
            bam = self.run_process(
                "upload-bam",
                {
                    "src": "bamclipper/output/TP53.primerclipped.bam",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )
            genome = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "bqsr/input/hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

        alignment_summary = self.run_process(
            "insert-size", {"bam": bam.id, "genome": genome.id}
        )

        self.assertFile(
            alignment_summary,
            "report",
            "hs_gatk_insert_size_metrics.txt",
            file_filter=filter_comment_lines,
        )

    @tag_process("wgs-metrics")
    def test_wgs_metrics(self):
        with self.preparation_stage():
            bam = self.run_process(
                "upload-bam",
                {
                    "src": "bamclipper/output/TP53.primerclipped.bam",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )
            genome = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "bqsr/input/hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

        wgs_metrics = self.run_process(
            "wgs-metrics",
            {"bam": bam.id, "genome": genome.id, "create_histogram": True},
        )

        self.assertFile(
            wgs_metrics,
            "report",
            "hs_gatk_wgs_metrics.txt",
            file_filter=filter_comment_lines,
        )

    @tag_process("rrbs-metrics")
    def test_rrbs_metrics(self):
        with self.preparation_stage():
            bam = self.run_process(
                "upload-bam",
                {
                    "src": "bamclipper/output/TP53.primerclipped.bam",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )
            genome = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "bqsr/input/hs_b37_chr17_upto_TP53.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

        rrbs_metrics = self.run_process(
            "rrbs-metrics",
            {
                "bam": bam.id,
                "genome": genome.id,
            },
        )

        self.assertFileExists(
            rrbs_metrics,
            "report",
        )

    @tag_process("merge-fastq-single", "merge-fastq-paired")
    def test_merge_fastq(self):
        self.collection.set_permission(Permission.VIEW, self.contributor)
        with self.preparation_stage():
            reads_single_1 = self.prepare_reads()
            reads_single_2 = self.prepare_reads()
            reads_paired_1 = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )
            reads_paired_2 = self.prepare_paired_reads()

            rel_type_group = RelationType.objects.get(name="group")

            replicate_group = Relation.objects.create(
                contributor=self.contributor,
                collection=self.collection,
                type=rel_type_group,
                category="Replicate",
            )

            treatment_group = Relation.objects.create(
                contributor=self.contributor,
                collection=self.collection,
                type=rel_type_group,
                category="Treatment",
            )

            RelationPartition.objects.create(
                relation=replicate_group,
                entity=reads_single_1.entity,
                label="single_sample",
            )
            RelationPartition.objects.create(
                relation=replicate_group,
                entity=reads_single_2.entity,
                label="single_sample",
            )
            RelationPartition.objects.create(
                relation=replicate_group,
                entity=reads_paired_1.entity,
                label="paired_sample",
            )
            RelationPartition.objects.create(
                relation=replicate_group,
                entity=reads_paired_2.entity,
                label="paired_sample",
            )

        inputs = {"reads": [reads_single_1.pk, reads_single_2.pk]}
        self.run_process("merge-fastq-single", inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        merged_1 = Data.objects.filter(process__slug="upload-fastq-single").last()
        self.assertFiles(
            merged_1,
            "fastq",
            [os.path.join("merge-fastq", "output", "reads_merged.fastq.gz")],
            compression="gzip",
        )

        RelationPartition.objects.create(
            relation=treatment_group,
            entity=reads_single_1.entity,
            label="treatment 1",
        )

        inputs = {"reads": [reads_single_1.pk, reads_single_2.pk]}
        merge_process = self.run_process("merge-fastq-single", inputs)

        relation_warn = (
            "Sample reads.fastq.gz has defined Treatment relation. "
            "Samples will only be merged based on Replicate relations."
        )
        self.assertEqual(
            merge_process.process_warning,
            [relation_warn],
        )

        inputs = {"reads": [reads_paired_1.pk, reads_paired_2.pk]}
        self.run_process("merge-fastq-paired", inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        merged_2 = Data.objects.filter(process__slug="upload-fastq-paired").last()
        self.assertFiles(
            merged_2,
            "fastq",
            [os.path.join("merge-fastq", "output", "fw_reads_merged.fastq.gz")],
            compression="gzip",
        )
        self.assertFiles(
            merged_2,
            "fastq2",
            [os.path.join("merge-fastq", "output", "rw_reads_merged.fastq.gz")],
            compression="gzip",
        )

    @tag_process("bedtools-bamtobed")
    def test_bedtools_bamtobed(self):
        with self.preparation_stage():
            species = "Dictyostelium discoideum"
            build = "dd-05-2009"

            bam = self.run_process(
                "upload-bam",
                {
                    "src": "reads.bam",
                    "species": species,
                    "build": build,
                },
            )

        bed = self.run_process(
            "bedtools-bamtobed",
            {
                "alignment": bam.id,
            },
        )

        self.assertFile(bed, "bedpe", Path("bam_processing") / "input" / "reads.bedpe")
        self.assertFields(bed, "species", species)
        self.assertFields(bed, "build", build)

    @tag_process("calculate-bigwig")
    def test_scale_bigwig(self):
        base = Path("bam_processing")
        inputs = base / "input"
        outputs = base / "output"

        with self.preparation_stage():
            species = "Dictyostelium discoideum"
            build = "dd-05-2009"

            bam = self.run_process(
                "upload-bam",
                {
                    "src": "reads.bam",
                    "species": species,
                    "build": build,
                },
            )

            bedpe = self.run_process(
                "upload-bedpe",
                {
                    "src": inputs / "reads.bedpe",
                    "species": species,
                    "build": build,
                },
            )

        bigwig = self.run_process("calculate-bigwig", {"alignment": bam.id})
        self.assertFile(bigwig, "bigwig", outputs / "reads.bigwig")
        self.assertFields(bigwig, "species", species)
        self.assertFields(bigwig, "build", build)

        scaled_bigwig = self.run_process(
            "calculate-bigwig",
            {
                "alignment": bam.id,
                "bedpe": bedpe.id,
                "scale": 10000,
            },
        )

        self.assertFile(scaled_bigwig, "bigwig", outputs / "reads.SInorm.bigwig")
        self.assertFields(scaled_bigwig, "species", species)
        self.assertFields(scaled_bigwig, "build", build)

    @tag_process("bamtofastq-paired")
    def test_bamtofastq_paired(self):
        base = Path("bamtofastq")
        inputs = base / "input"
        outputs = base / "output"
        with self.preparation_stage():
            bam = self.run_process(
                "upload-bam",
                {
                    "src": inputs / "alignment.bam",
                    "species": "Homo sapiens",
                    "build": "hg38",
                },
            )

        reads = self.run_process(
            "bamtofastq-paired",
            {"bam": bam.id},
        )

        self.assertFiles(
            reads,
            "fastq",
            [outputs / "output_reads_mate1.fastq.gz"],
            compression="gzip",
        )

        self.assertFiles(
            reads,
            "fastq2",
            [outputs / "output_reads_mate2.fastq.gz"],
            compression="gzip",
        )

    @tag_process("mutations-table")
    def test_report_variants(self):
        input_folder = Path("report_variants") / "input"
        output_folder = Path("report_variants") / "output"
        with self.preparation_stage():
            variants_vcf = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "KRAS_mutation.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            variants_clinvar = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "filtered_variants_clinvar.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            no_mutation = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "no_mutation.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            no_input = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "empty.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            ensembl = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "ensembl.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            snpeff_ensembl = self.run_process(
                "snpeff-single",
                {
                    "variants": ensembl.id,
                },
            )
            snpeff = self.run_process(
                "snpeff-single",
                {
                    "variants": variants_vcf.id,
                },
            )
            snpeff_clinvar = self.run_process(
                "snpeff-single", {"variants": variants_clinvar.id}
            )
            geneset = self.run_process(
                "create-geneset",
                {
                    "genes": ["ENSG00000133703"],
                    "species": "Homo sapiens",
                    "source": "ENSEMBL",
                },
            )
            snpeff_nomutation = self.run_process(
                "snpeff-single",
                {
                    "variants": no_mutation.id,
                },
            )
            snpeff_noinput = self.run_process(
                "snpeff-single",
                {
                    "variants": no_input.id,
                },
            )
            bam = self.run_process(
                "upload-bam",
                {
                    "src": input_folder / "chr12_KRAS_sorted.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            ensembl_bam = self.run_process(
                "upload-bam",
                {
                    "src": input_folder / "ensembl.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

        snpeff.entity = bam.entity
        snpeff_nomutation.entity = bam.entity
        snpeff_noinput.entity = bam.entity
        snpeff_clinvar.entity = bam.entity
        snpeff_ensembl.entity = ensembl_bam.entity
        snpeff_ensembl.save()

        snpeff.save()
        snpeff_nomutation.save()
        snpeff_noinput.save()
        snpeff_clinvar.save()

        report = self.run_process(
            "mutations-table",
            {
                "variants": snpeff.id,
                "bam": bam.id,
                "mutations": ["KRAS: Gly12, Gly13"],
                "advanced": {"gf_fields": ["DP"]},
            },
        )
        self.assertFile(report, "tsv", output_folder / "mutations_KRAS.tsv")

        experiment = VariantExperiment.objects.filter().last()
        self.assertEqual(
            experiment.variant_data_source, "RNA-Seq Variant Calling Pipeline"
        )

        variant = Variant.objects.get(
            species="Homo sapiens",
            genome_assembly="GRCh38",
            chromosome="chr12",
            position=25245347,
            reference="C",
            alternative="T",
        )

        variant_call = VariantCall.objects.filter(variant=variant)[0]
        self.assertEqual(variant_call.sample_id, report.entity.id)
        self.assertEqual(variant_call.data_id, report.id)
        self.assertEqual(variant_call.quality, 974.64)
        self.assertEqual(variant_call.depth, 41)
        self.assertEqual(variant_call.depth_norm_quality, 23.77)
        self.assertEqual(variant_call.filter, "PASS")
        self.assertEqual(variant_call.unfiltered_allele_depth, None)
        self.assertEqual(variant_call.genotype, None)
        self.assertEqual(variant_call.genotype_quality, None)

        report = self.run_process(
            "mutations-table",
            {
                "variants": snpeff_nomutation.id,
                "bam": bam.id,
                "mutations": ["KRAS"],
                "vcf_fields": [
                    "CHROM",
                    "POS",
                    "ID",
                    "QUAL",
                    "REF",
                    "ALT",
                    "FILTER",
                    "ANN",
                ],
                "advanced": {"gf_fields": ["GT", "AD", "DP"]},
            },
        )
        self.assertFile(report, "tsv", output_folder / "no_mutations.tsv")

        report = self.run_process(
            "mutations-table",
            {
                "variants": snpeff_noinput.id,
                "bam": bam.id,
                "mutations": ["KRAS: Gly12, Gly13"],
            },
        )
        self.assertEqual(
            report.process_warning,
            [
                "There are no variants in the input VCF file.",
            ],
        )

        report = self.run_process(
            "mutations-table",
            {
                "variants": snpeff.id,
                "vcf_fields": ["CHROM", "POS", "ID", "QUAL", "REF", "ALT", "ANN"],
                "bam": bam.id,
                "mutations": ["KRAS: Ala11"],
            },
        )
        self.assertEqual(
            report.process_warning,
            ["No variants present for the input set of mutations."],
        )

        report = self.run_process(
            "mutations-table",
            {
                "variants": snpeff_nomutation.id,
                "bam": bam.id,
                "mutations": ["KRAS: Gly12,Gly13, Gln61, Ar3"],
            },
            Data.STATUS_ERROR,
        )
        self.assertEqual(
            report.process_error,
            [
                "The input amino acid Ar3 is in the wrong format "
                "or is not among the 20 standard amino acids."
            ],
        )

        report = self.run_process(
            "mutations-table",
            {
                "variants": snpeff_clinvar.id,
                "bam": bam.id,
                "geneset": geneset.id,
                "vcf_fields": [
                    "CHROM",
                    "POS",
                    "ID",
                    "QUAL",
                    "REF",
                    "ALT",
                    "DP",
                    "ANN",
                    "CLNDN",
                    "CLNSIG",
                ],
            },
        )
        self.assertFile(
            report,
            "tsv",
            output_folder / "mutations_clinvar.tsv",
        )

        ensembl_mutations = self.run_process(
            "mutations-table",
            {
                "variants": snpeff_ensembl.id,
                "bam": ensembl_bam.id,
                "mutations": ["PTPRD"],
                "vcf_fields": [
                    "CHROM",
                    "POS",
                    "ID",
                    "QUAL",
                    "REF",
                    "ALT",
                    "FILTER",
                    "ANN",
                    "CLNDN",
                    "CLNSIG",
                    "QD",
                ],
                "advanced": {"gf_fields": ["GT", "GQ", "AD", "DP", "FT"]},
            },
        )
        self.assertFile(
            ensembl_mutations,
            "tsv",
            output_folder / "ensembl_variants.tsv",
        )

        variant_calls = VariantCall.objects.filter(
            sample_id=ensembl_mutations.entity.id
        )
        self.assertEqual(len(variant_calls), 3)

        variant_call_1 = variant_calls[0]
        self.assertEqual(variant_call_1.depth, 2)
        self.assertEqual(variant_call_1.quality, 47.32)
        self.assertEqual(variant_call_1.genotype, "1/1")
        self.assertEqual(variant_call_1.genotype_quality, 6)
        self.assertEqual(variant_call_1.unfiltered_allele_depth, 2)
        self.assertEqual(variant_call_1.depth_norm_quality, 23.66)
        self.assertEqual(variant_call_1.filter, "DP")
