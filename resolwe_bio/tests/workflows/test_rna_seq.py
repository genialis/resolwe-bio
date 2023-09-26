from pathlib import Path

from django.test import override_settings

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_docker_executor, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class RNASeqWorkflowTestCase(KBBioProcessTestCase):
    @tag_process("workflow-rnaseq-cuffquant")
    def test_cuffquant_workflow(self):
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn="genome.fasta.gz",
                species="Dictyostelium discoideum",
                build="dd-05-2009",
            )
            hisat2_index = self.run_process("hisat2-index", {"ref_seq": ref_seq.id})
            reads = self.prepare_reads()
            annotation = self.prepare_annotation_gff()

        self.run_process(
            "workflow-rnaseq-cuffquant",
            {
                "reads": reads.id,
                "genome": hisat2_index.id,
                "annotation": annotation.id,
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

    @with_resolwe_host
    @tag_process("workflow-quantseq")
    def test_bbduk_star_fc_quant_workflow(self):
        output_folder = Path("quantseq") / "output"
        with self.preparation_stage():
            inputs = {"src": ["hs_single bbduk_star_htseq_reads_single.fastq.gz"]}
            reads = self.run_processor("upload-fastq-single", inputs)

            paired_reads = self.prepare_paired_reads(
                ["hs_paired_R1 workflow_bbduk_star_htseq.fastq.gz"],
                ["hs_paired_R2 workflow_bbduk_star_htseq.fastq.gz"],
            )

            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "hs genome.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )
            adapters = self.prepare_ref_seq()

            inputs = {
                "src": "hs annotation.gtf.gz",
                "source": "ENSEMBL",
                "species": "Homo sapiens",
                "build": "ens_90",
            }
            annotation = self.run_process("upload-gtf", inputs)

            inputs = {"annotation": annotation.id, "ref_seq": star_index_fasta.id}
            star_index = self.run_process("alignment-star-index", inputs)

            rrna_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "Homo_sapiens_rRNA.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "rRNA",
                },
            )
            rrna_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": rrna_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

            globin_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "Homo_sapiens_globin_reference.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "globin",
                },
            )
            globin_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": globin_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

        self.run_process(
            "workflow-quantseq",
            {
                "trimming_tool": "bbduk",
                "reads": reads.id,
                "adapters": [adapters.id],
                "genome": star_index.id,
                "annotation": annotation.id,
                "assay_type": "forward",
                "rrna_reference": rrna_star_index.id,
                "globin_reference": globin_star_index.id,
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        feature_counts = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            feature_counts,
            "rc",
            output_folder / "workflow_bbduk_star_fc_quant_single_rc.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            feature_counts,
            "exp",
            output_folder / "workflow_bbduk_star_fc_quant_single_cpm.tab.gz",
            compression="gzip",
        )
        self.assertFields(feature_counts, "exp_type", "CPM")
        self.assertFields(feature_counts, "source", "ENSEMBL")
        self.assertFields(feature_counts, "species", "Homo sapiens")

        self.run_process(
            "workflow-quantseq",
            {
                "trimming_tool": "bbduk",
                "reads": paired_reads.id,
                "adapters": [adapters.id],
                "genome": star_index.id,
                "annotation": annotation.id,
                "assay_type": "reverse",
                "rrna_reference": rrna_star_index.id,
                "globin_reference": globin_star_index.id,
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        fc_paired = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            fc_paired,
            "rc",
            output_folder / "workflow_bbduk_star_fc_quant_paired_rc.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            fc_paired,
            "exp",
            output_folder / "workflow_bbduk_star_fc_quant_paired_cpm.tab.g",
            compression="gzip",
        )
        self.assertFields(fc_paired, "exp_type", "CPM")
        self.assertFields(fc_paired, "source", "ENSEMBL")
        self.assertFields(fc_paired, "species", "Homo sapiens")

    @with_resolwe_host
    @tag_process("workflow-quantseq")
    def test_cutadapt_star_fc_quant_workflow(self):
        output_folder = Path("quantseq") / "output"
        with self.preparation_stage():
            reads = self.run_processor(
                "upload-fastq-single",
                {"src": ["hs_single bbduk_star_htseq_reads_single.fastq.gz"]},
            )

            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "hs genome.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            annotation = self.run_process(
                "upload-gtf",
                {
                    "src": "hs annotation.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            star_index = self.run_process(
                "alignment-star-index",
                {"annotation": annotation.id, "ref_seq": star_index_fasta.id},
            )

            rrna_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "Homo_sapiens_rRNA.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "rRNA",
                },
            )

            rrna_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": rrna_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

            globin_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "Homo_sapiens_globin_reference.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "globin",
                },
            )

            globin_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": globin_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

        self.run_process(
            "workflow-quantseq",
            {
                "trimming_tool": "cutadapt",
                "reads": reads.id,
                "genome": star_index.id,
                "annotation": annotation.id,
                "rrna_reference": rrna_star_index.id,
                "globin_reference": globin_star_index.id,
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        feature_counts = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            feature_counts,
            "rc",
            output_folder / "workflow_cutadapt_star_fc_quant_single_rc.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            feature_counts,
            "exp",
            output_folder / "workflow_cutadapt_star_fc_quant_single_cpm.tab.gz",
            compression="gzip",
        )
        self.assertFields(feature_counts, "exp_type", "CPM")
        self.assertFields(feature_counts, "source", "ENSEMBL")
        self.assertFields(feature_counts, "species", "Homo sapiens")

        self.run_process(
            "workflow-quantseq",
            {
                "trimming_tool": "cutadapt",
                "reads": reads.id,
                "genome": star_index.id,
                "annotation": annotation.id,
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        feature_counts = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            feature_counts,
            "rc",
            output_folder / "workflow_cutadapt_star_fc_quant_single_rc.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            feature_counts,
            "exp",
            output_folder / "workflow_cutadapt_star_fc_quant_single_cpm.tab.gz",
            compression="gzip",
        )
        self.assertFields(feature_counts, "exp_type", "CPM")
        self.assertFields(feature_counts, "source", "ENSEMBL")
        self.assertFields(feature_counts, "species", "Homo sapiens")

    @with_resolwe_host
    @with_docker_executor
    @override_settings(FLOW_PROCESS_MAX_CORES=4)
    @tag_process(
        "workflow-bbduk-star-featurecounts-qc",
    )
    def test_bbduk_star_featurecounts_workflow(self):
        input_folder = Path("test_star") / "input"
        output_folder = Path("test_featurecounts") / "outputs"
        with self.preparation_stage():
            reads = self.prepare_reads(["chr1_single.fastq.gz"])
            paired_reads = self.prepare_paired_reads(
                ["chr1_paired_R1.fastq.gz"], ["chr1_paired_R2.fastq.gz"]
            )
            paired_lanes = self.prepare_paired_reads(
                mate1=[
                    input_folder / "chr1_paired_R1_mate1.fastq.gz",
                    input_folder / "chr1_paired_R1_mate2.fastq.gz",
                ],
                mate2=[
                    input_folder / "chr1_paired_R2_mate1.fastq.gz",
                    input_folder / "chr1_paired_R2_mate2.fastq.gz",
                ],
            )
            annotation = self.prepare_annotation(
                fn="chr1_region.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="ens90",
            )
            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "chr1_region.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens90",
                },
            )
            inputs = {
                "annotation": annotation.id,
                "ref_seq": star_index_fasta.id,
            }
            star_index = self.run_process("alignment-star-index", inputs)
            adapters = self.prepare_ref_seq()

            rrna_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "Homo_sapiens_rRNA.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "rRNA",
                },
            )
            rrna_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": rrna_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

            globin_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "Homo_sapiens_globin_reference.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "globin",
                },
            )
            globin_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": globin_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

        inputs = {
            "reads": reads.id,
            "genome": star_index.id,
            "annotation": annotation.id,
            "rrna_reference": rrna_star_index.id,
            "globin_reference": globin_star_index.id,
            "preprocessing": {
                "adapters": [adapters.id],
                "custom_adapter_sequences": ["ACTGACTGACTG", "AAACCCTTT"],
            },
        }

        self.run_process("workflow-bbduk-star-featurecounts-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        feature_counts = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            feature_counts, "rc", "feature_counts_rc_single.tab.gz", compression="gzip"
        )
        self.assertEqual(feature_counts.name, "Quantified (chr1_single.fastq.gz)")

        downsampled = Data.objects.filter(process__slug="alignment-star").last()
        self.assertFields(downsampled, "build", "ens90")
        self.assertEqual(downsampled.name, "Aligned subset (chr1_single.fastq.gz)")

        multiqc = Data.objects.filter(process__slug="multiqc").last()
        self.assertFileExists(multiqc, "report")

        inputs["reads"] = paired_reads.id
        self.run_process("workflow-bbduk-star-featurecounts-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        feature_counts = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            feature_counts,
            "rc",
            "feature_counts_rc_paired_2.tab.gz",
            compression="gzip",
        )
        self.assertEqual(feature_counts.name, "Quantified (chr1_paired_R1.fastq.gz)")

        downsampled = Data.objects.filter(process__slug="alignment-star").last()
        self.assertFields(downsampled, "build", "ens90")
        self.assertEqual(downsampled.name, "Aligned subset (chr1_paired_R1.fastq.gz)")

        multiqc = Data.objects.filter(process__slug="multiqc").last()
        self.assertFileExists(multiqc, "report")

        # test the pipeline without the adapter sequences specified
        del inputs["preprocessing"]["adapters"]
        del inputs["preprocessing"]["custom_adapter_sequences"]
        self.run_process("workflow-bbduk-star-featurecounts-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        feature_counts = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            feature_counts,
            "rc",
            "feature_counts_rc_paired_wo_adapter_trim_2.tab.gz",
            compression="gzip",
        )
        self.assertEqual(feature_counts.name, "Quantified (chr1_paired_R1.fastq.gz)")

        inputs = {
            "reads": paired_lanes.id,
            "genome": star_index.id,
            "annotation": annotation.id,
            "rrna_reference": rrna_star_index.id,
            "globin_reference": globin_star_index.id,
            "preprocessing": {
                "adapters": [adapters.id],
                "custom_adapter_sequences": ["ACTGACTGACTG", "AAACCCTTT"],
            },
        }
        self.run_process("workflow-bbduk-star-featurecounts-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        feature_counts = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            feature_counts,
            "exp_set",
            "chr1_paired_R1_workflow_bbduk_star_htseq_preprocessed_expressions.txt.gz",
            compression="gzip",
        )
        self.assertFile(
            feature_counts,
            "per_lane_rc",
            output_folder / "per_lane_rc.txt.gz",
            compression="gzip",
        )

    @with_resolwe_host
    @tag_process("workflow-corall-single", "workflow-corall-paired")
    def test_corall(self):
        with self.preparation_stage():
            reads = self.prepare_reads(["./corall/input/corall_single.fastq.gz"])

            reads_paired = self.prepare_paired_reads(
                mate1=["./corall/input/corall_mate1.fastq.gz"],
                mate2=["./corall/input/corall_mate2.fastq.gz"],
            )

            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "./corall/input/hs_genome_chr2_1_45000.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            annotation = self.run_process(
                "upload-gtf",
                {
                    "src": "./corall/input/hs_annotation_chr2_1_45000.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            star_index = self.run_process(
                "alignment-star-index",
                {
                    "annotation": annotation.id,
                    "ref_seq": star_index_fasta.id,
                },
            )

            rrna_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "./corall/input/Homo_sapiens_rRNA.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "rRNA",
                },
            )

            rrna_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": rrna_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

            globin_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "./corall/input/Homo_sapiens_globin_reference.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "globin",
                },
            )

            globin_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": globin_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

        self.run_process(
            "workflow-corall-single",
            {
                "reads": reads.id,
                "star_index": star_index.id,
                "annotation": annotation.id,
                "rrna_reference": rrna_star_index.id,
                "globin_reference": globin_star_index.id,
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        exp = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            exp,
            "exp_set",
            "./corall/output/corall_workfow_expression_single.txt.gz",
            compression="gzip",
        )
        self.assertFields(exp, "exp_type", "TPM")
        self.assertFields(exp, "source", "ENSEMBL")
        self.assertFields(exp, "species", "Homo sapiens")

        self.run_process(
            "workflow-corall-paired",
            {
                "reads": reads_paired.id,
                "star_index": star_index.id,
                "annotation": annotation.id,
                "rrna_reference": rrna_star_index.id,
                "globin_reference": globin_star_index.id,
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        exp = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            exp,
            "exp_set",
            "./corall/output/corall_workfow_expression_paired.txt.gz",
            compression="gzip",
        )
        self.assertFields(exp, "exp_type", "TPM")
        self.assertFields(exp, "source", "ENSEMBL")
        self.assertFields(exp, "species", "Homo sapiens")

    @with_resolwe_host
    @tag_process("workflow-bbduk-salmon-qc")
    def test_salmon_workflow(self):
        with self.preparation_stage():
            reads = self.prepare_reads(["chr1_single.fastq.gz"])
            paired_reads = self.prepare_paired_reads(
                ["chr1_paired_R1.fastq.gz"], ["chr1_paired_R2.fastq.gz"]
            )
            annotation = self.prepare_annotation(
                fn="chr1_region_salmon.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="ens92",
            )
            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "chr1_region.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens92",
                },
            )
            star_index = self.run_process(
                "alignment-star-index",
                {
                    "annotation": annotation.id,
                    "ref_seq": star_index_fasta.id,
                },
            )
            cdna = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "chr1_region.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens92",
                },
            )
            salmon_index = self.run_process(
                "salmon-index",
                {
                    "nucl": cdna.id,
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens92",
                },
            )
            adapters = self.prepare_ref_seq()
            rrna_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "salmon_workflow/input/Homo_sapiens_rRNA.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "rRNA",
                },
            )
            rrna_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": rrna_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )
            globin_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "salmon_workflow/input/Homo_sapiens_globin_reference.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "globin",
                },
            )
            globin_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": globin_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

        inputs = {
            "reads": reads.id,
            "genome": star_index.id,
            "salmon_index": salmon_index.id,
            "annotation": annotation.id,
            "rrna_reference": rrna_star_index.id,
            "globin_reference": globin_star_index.id,
            "preprocessing": {
                "adapters": [adapters.id],
            },
            "quantification": {
                "min_assigned_frag": 1,
                "gc_bias": False,
            },
        }

        self.run_process("workflow-bbduk-salmon-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        salmon_single_end = Data.objects.get(process__slug="salmon-quant")
        self.assertFile(
            salmon_single_end,
            "exp_set",
            "salmon_workflow/output/single_end_exp.txt.gz",
            compression="gzip",
        )
        self.assertFields(salmon_single_end, "exp_type", "TPM")
        self.assertFields(salmon_single_end, "source", "ENSEMBL")

        inputs["reads"] = paired_reads.id
        inputs["quantification"]["gc_bias"] = True
        inputs["quantification"]["num_bootstraps"] = 5
        self.run_process("workflow-bbduk-salmon-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        salmon_paired_end = Data.objects.filter(process__slug="salmon-quant").last()
        self.assertFile(
            salmon_paired_end,
            "exp_set",
            "salmon_workflow/output/paired_end_exp.txt.gz",
            compression="gzip",
        )
        self.assertFields(salmon_paired_end, "exp_type", "TPM")
        self.assertFields(salmon_paired_end, "source", "ENSEMBL")
        self.assertFileExists(salmon_paired_end, "variance")


class RNASeqVCWorkflowTestCase(KBBioProcessTestCase):
    @tag_process("workflow-rnaseq-variantcalling")
    def test_rnaseq_variantcalling(self):
        input_folder = Path("rnaseq_variantcalling") / "input"
        output_folder = Path("rnaseq_variantcalling") / "output"
        with self.preparation_stage():
            reads = self.run_process(
                "upload-fastq-single",
                {"src": [input_folder / "chr1_19000_R1.fastq.gz"]},
            )
            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": input_folder / "chr1_19000.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            adapters = self.prepare_ref_seq()
            star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": ref_seq.id,
                    "source": "ENSEMBL",
                },
            )
            alignment = self.run_process(
                "alignment-star",
                {
                    "reads": reads.id,
                    "genome": star_index.id,
                    "two_pass_mapping": {"two_pass_mode": True},
                    "output_options": {"out_unmapped": True},
                },
            )
            alignment_onepass = self.run_process(
                "alignment-star",
                {
                    "reads": reads.id,
                    "genome": star_index.id,
                },
            )
            dbsnp = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "dbsnp-hg38.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            intervals = self.run_process(
                "upload-bed",
                {
                    "src": input_folder / "hg38.intervals.bed",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            geneset = self.run_process(
                "create-geneset",
                {
                    "genes": ["ENSG00000223972"],
                    "species": "Homo sapiens",
                    "source": "ENSEMBL",
                },
            )

        input_workflow = {
            "bam": alignment.id,
            "ref_seq": ref_seq.id,
            "dbsnp": dbsnp.id,
            "mutations": ["DDX11L1"],
            "variant_filtration": {"mask": dbsnp.id, "mask_name": "DB"},
        }

        self.run_process(
            process_slug="workflow-rnaseq-variantcalling", input_=input_workflow
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        mutations = Data.objects.last()
        self.assertFile(mutations, "tsv", output_folder / "mutations_bam.tsv")

        input_workflow = {
            "bam": alignment_onepass.id,
            "ref_seq": ref_seq.id,
            "dbsnp": dbsnp.id,
            "mutations": ["DDX11L1"],
            "advanced": {"multiqc": True},
        }

        warning_onepass = self.run_process(
            process_slug="workflow-rnaseq-variantcalling", input_=input_workflow
        )

        warning_msg = [
            "Two-pass mode was not used in alignment with STAR. It is "
            "highly recommended that you use two-pass mode for RNA-seq "
            "variant calling.",
            "It is recommended that you use parameter '--outSAMunmapped Within' "
            "in STAR alignment.",
        ]
        self.assertEqual(warning_onepass.process_warning, warning_msg)

        multiqc = Data.objects.filter(process__slug="multiqc").last()
        self.assertFileExists(multiqc, "report")

        input_workflow = {
            "reads": reads.id,
            "bbduk": {
                "adapters": [adapters.id],
            },
            "ref_seq": ref_seq.id,
            "genome": star_index.id,
            "dbsnp": dbsnp.id,
            "mutations": ["DDX11L1"],
        }

        self.run_process(
            process_slug="workflow-rnaseq-variantcalling", input_=input_workflow
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        mutations = Data.objects.filter(process__slug="mutations-table").last()
        self.assertFile(mutations, "tsv", output_folder / "mutations.tsv")

        # Test for workflow without reads preprocessing
        input_workflow = {
            "reads": reads.id,
            "preprocessing": False,
            "ref_seq": ref_seq.id,
            "genome": star_index.id,
            "dbsnp": dbsnp.id,
            "intervals": intervals.id,
            "geneset": geneset.id,
            "clinvar": dbsnp.id,
            "haplotype_caller": {"interval_padding": 0},
        }
        self.run_process(
            process_slug="workflow-rnaseq-variantcalling", input_=input_workflow
        )
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        mutations = Data.objects.filter(process__slug="mutations-table").last()
        self.assertFile(mutations, "tsv", output_folder / "mutations_geneset.tsv")


class STARRNASeqWorkflowTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @with_docker_executor
    @override_settings(FLOW_PROCESS_MAX_CORES=4)
    @tag_process("workflow-bbduk-star-qc")
    def test_bbduk_star_workflow(self):
        output_folder = Path("test_star") / "output"
        input_folder = Path("test_star") / "input"
        with self.preparation_stage():
            reads = self.prepare_reads(["chr1_single.fastq.gz"])
            paired_reads = self.prepare_paired_reads(
                ["chr1_paired_R1.fastq.gz"], ["chr1_paired_R2.fastq.gz"]
            )
            paired_lanes = self.prepare_paired_reads(
                mate1=[
                    input_folder / "chr1_paired_R1_mate1.fastq.gz",
                    input_folder / "chr1_paired_R1_mate2.fastq.gz",
                ],
                mate2=[
                    input_folder / "chr1_paired_R2_mate1.fastq.gz",
                    input_folder / "chr1_paired_R2_mate2.fastq.gz",
                ],
            )
            annotation = self.prepare_annotation(
                fn="chr1_region.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="ens90",
            )
            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "chr1_region.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens90",
                },
            )
            cdna = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "chr1_region.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "ens92",
                },
            )
            salmon_index = self.run_process(
                "salmon-index",
                {
                    "nucl": cdna.id,
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens92",
                },
            )
            inputs = {
                "annotation": annotation.id,
                "ref_seq": star_index_fasta.id,
            }
            star_index = self.run_process("alignment-star-index", inputs)
            adapters = self.prepare_ref_seq()

            rrna_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "Homo_sapiens_rRNA.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "rRNA",
                },
            )
            rrna_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": rrna_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

            globin_reference = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "Homo_sapiens_globin_reference.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "globin",
                },
            )
            globin_star_index = self.run_process(
                "alignment-star-index",
                {
                    "ref_seq": globin_reference.id,
                    "source": "NCBI",
                    "advanced": {
                        "genome_sa_string_len": 2,
                    },
                },
            )

        inputs = {
            "reads": reads.id,
            "genome": star_index.id,
            "annotation": annotation.id,
            "rrna_reference": rrna_star_index.id,
            "globin_reference": globin_star_index.id,
            "preprocessing": {
                "adapters": [adapters.id],
                "custom_adapter_sequences": ["ACTGACTGACTG", "AAACCCTTT"],
            },
        }

        self.run_process("workflow-bbduk-star-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        star_quant = Data.objects.filter(process__slug="star-quantification").last()
        self.assertFile(
            star_quant, "rc", "feature_counts_rc_single.tab.gz", compression="gzip"
        )
        self.assertEqual(star_quant.name, "Quantified (chr1_single.fastq.gz)")

        multiqc = Data.objects.filter(process__slug="multiqc").last()
        self.assertFileExists(multiqc, "report")

        inputs["reads"] = paired_reads.id
        self.run_process("workflow-bbduk-star-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        star_quant = Data.objects.filter(process__slug="star-quantification").last()
        self.assertFile(
            star_quant, "rc", "feature_counts_rc_paired.tab.gz", compression="gzip"
        )
        self.assertEqual(star_quant.name, "Quantified (chr1_paired_R1.fastq.gz)")

        multiqc = Data.objects.filter(process__slug="multiqc").last()
        self.assertFileExists(multiqc, "report")

        # test the pipeline without the adapter sequences specified
        del inputs["preprocessing"]["adapters"]
        del inputs["preprocessing"]["custom_adapter_sequences"]
        self.run_process("workflow-bbduk-star-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        star_quant = Data.objects.filter(process__slug="star-quantification").last()
        self.assertFile(
            star_quant,
            "rc",
            "feature_counts_rc_paired_wo_adapter_trim.tab.gz",
            compression="gzip",
        )
        self.assertEqual(star_quant.name, "Quantified (chr1_paired_R1.fastq.gz)")

        inputs = {
            "reads": paired_lanes.id,
            "genome": star_index.id,
            "annotation": annotation.id,
            "rrna_reference": rrna_star_index.id,
            "globin_reference": globin_star_index.id,
            "preprocessing": {
                "adapters": [adapters.id],
                "custom_adapter_sequences": ["ACTGACTGACTG", "AAACCCTTT"],
            },
        }
        self.run_process("workflow-bbduk-star-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        star_quant = Data.objects.filter(process__slug="star-quantification").last()
        self.assertFile(
            star_quant,
            "exp_set",
            output_folder / "workflow_expressions.txt.gz",
            compression="gzip",
        )

        inputs["assay_type"] = "non_specific"
        inputs["cdna_index"] = salmon_index.id
        inputs.update({"alignment": {"two_pass_mapping": {"two_pass_mode": False}}})
        self.run_process("workflow-bbduk-star-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        star_quant = Data.objects.filter(process__slug="star-quantification").last()
        self.assertFile(
            star_quant,
            "exp_set",
            output_folder / "workflow_expressions_2.txt.gz",
            compression="gzip",
        )
        multiqc = Data.objects.filter(process__slug="multiqc").last()
        self.assertFileExists(multiqc, "report")
