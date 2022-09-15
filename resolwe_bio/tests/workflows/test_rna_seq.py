from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.filter import filter_vcf_variable
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
    @tag_process(
        "workflow-bbduk-star-featurecounts-qc",
    )
    def test_bbduk_star_featurecounts_workflow(self):
        with self.preparation_stage():
            reads = self.prepare_reads(["hs sim_reads_single.fastq.gz"])
            paired_reads = self.prepare_paired_reads(
                ["hs sim_reads1.fastq.gz"], ["hs sim_reads2.fastq.gz"]
            )
            annotation = self.prepare_annotation(
                fn="hs annotation.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="ens90",
            )
            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "hs genome.fasta.gz",
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
        self.assertEqual(
            feature_counts.name, "Quantified (hs sim_reads_single.fastq.gz)"
        )

        globin = Data.objects.filter(process__slug="alignment-star").last()
        self.assertFields(globin, "build", "globin")
        self.assertEqual(globin.name, "Globin aligned (hs sim_reads_single.fastq.gz)")

        multiqc = Data.objects.filter(process__slug="multiqc").last()
        self.assertFileExists(multiqc, "report")

        inputs["reads"] = paired_reads.id
        self.run_process("workflow-bbduk-star-featurecounts-qc", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        feature_counts = Data.objects.filter(process__slug="feature_counts").last()
        self.assertFile(
            feature_counts, "rc", "feature_counts_rc_paired.tab.gz", compression="gzip"
        )
        self.assertEqual(feature_counts.name, "Quantified (hs sim_reads1.fastq.gz)")

        globin = Data.objects.filter(process__slug="alignment-star").last()
        self.assertFields(globin, "build", "globin")
        self.assertEqual(globin.name, "Globin aligned (hs sim_reads1.fastq.gz)")

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
            "feature_counts_rc_paired_wo_adapter_trim.tab.gz",
            compression="gzip",
        )
        self.assertEqual(feature_counts.name, "Quantified (hs sim_reads1.fastq.gz)")

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
            reads = self.prepare_reads(
                ["salmon_workflow/input/hs sim_reads_single.fastq.gz"]
            )
            paired_reads = self.prepare_paired_reads(
                ["salmon_workflow/input/hs sim_reads1.fastq.gz"],
                ["salmon_workflow/input/hs sim_reads2.fastq.gz"],
            )
            annotation = self.prepare_annotation(
                fn="salmon_workflow/input/hs annotation.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="ens92",
            )
            star_index_fasta = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": "salmon_workflow/input/hs genome.fasta.gz",
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
                    "src": "salmon_workflow/input/hs cdna.fasta.gz",
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

    @tag_process(
        "workflow-rnaseq-variantcalling", "workflow-rnaseq-variantcalling-beta"
    )
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
                    "build": "custom_build",
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

            dbsnp = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "dbsnp-hg38.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            intervals = self.run_process(
                "upload-bed",
                {
                    "src": input_folder / "hg38.intervals.bed",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

        input_workflow = {
            "reads": reads.id,
            "bbduk": {
                "adapters": [adapters.id],
            },
            "ref_seq": ref_seq.id,
            "genome": star_index.id,
            "dbsnp": dbsnp.id,
            "exclude_filtered": True,
            "select_variants": {
                "select_type": ["SNP", "INDEL"],
            },
        }

        self.run_process(
            process_slug="workflow-rnaseq-variantcalling", input_=input_workflow
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        variants = Data.objects.filter(
            process__slug="gatk-select-variants-single"
        ).last()
        self.assertFile(
            variants,
            "vcf",
            output_folder / "selected_variants.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

        self.run_process(
            process_slug="workflow-rnaseq-variantcalling-beta", input_=input_workflow
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        variants = Data.objects.filter(
            process__slug="gatk-select-variants-single"
        ).last()
        self.assertFile(
            variants,
            "vcf",
            output_folder / "selected_variants.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

        # Test for workflow without reads preprocessing
        input_workflow["preprocessing"] = False
        input_workflow["intervals"] = intervals.id
        self.run_process(
            process_slug="workflow-rnaseq-variantcalling", input_=input_workflow
        )
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        variants = Data.objects.filter(
            process__slug="gatk-select-variants-single"
        ).last()
        self.assertFile(
            variants,
            "vcf",
            output_folder / "selected_variants_interval.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )
