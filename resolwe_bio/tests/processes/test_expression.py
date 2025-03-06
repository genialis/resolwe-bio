import os
from pathlib import Path

from resolwe.flow.models import Collection, Data, Relation
from resolwe.flow.models.entity import RelationPartition, RelationType
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.expression_filters.relation import replicate_groups
from resolwe_bio.utils.test import KBBioProcessTestCase


class ExpressionProcessorTestCase(KBBioProcessTestCase):
    fixtures = ["relationtypes.yaml"]

    @tag_process("cufflinks")
    def test_cufflinks(self):
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn="genome.fasta.gz",
                species="Dictyostelium discoideum",
                build="dd-05-2009",
            )
            hisat2_index = self.run_process("hisat2-index", {"ref_seq": ref_seq.id})
            reads = self.prepare_reads()
            annotation_gtf = self.prepare_annotation("annotation dicty.gff.gz")

            aligned_reads = self.run_process(
                "alignment-hisat2",
                {
                    "genome": hisat2_index.pk,
                    "reads": reads.pk,
                    "spliced_alignments": {"cufflinks": True},
                },
            )

        inputs = {
            "alignment": aligned_reads.pk,
            "annotation": annotation_gtf.pk,
            "genome": ref_seq.pk,
        }
        cuff_exp = self.run_process("cufflinks", inputs)
        self.assertFile(cuff_exp, "transcripts", "cufflinks_transcripts.gtf", sort=True)
        self.assertFields(cuff_exp, "species", "Dictyostelium discoideum")
        self.assertFields(cuff_exp, "build", "dd-05-2009")
        self.assertFields(cuff_exp, "source", "DICTYBASE")

    @tag_process("cuffquant")
    def test_cuffquant(self):
        with self.preparation_stage():
            inputs = {
                "src": "cuffquant_mapping.bam",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            bam = self.run_process("upload-bam", inputs)

            annotation = self.prepare_annotation(
                fn="hg19_chr20_small.gtf.gz",
                source="UCSC",
                species="Homo sapiens",
                build="hg19",
            )

        inputs = {"alignment": bam.id, "annotation": annotation.id}
        cuffquant = self.run_process("cuffquant", inputs)
        self.assertFields(cuffquant, "species", "Homo sapiens")
        self.assertFields(cuffquant, "build", "hg19")
        self.assertFields(cuffquant, "source", "UCSC")

    @with_resolwe_host
    @tag_process("cuffnorm")
    def test_cuffnorm(self):
        with self.preparation_stage():
            collection = Collection.objects.create(
                name="Test collection", contributor=self.contributor
            )

            rel_type_group = RelationType.objects.get(name="group")

            replicate_group = Relation.objects.create(
                contributor=self.contributor,
                collection=collection,
                type=rel_type_group,
                category="Replicate",
            )

            inputs = {
                "src": "cuffquant 1.cxb",
                "source": "UCSC",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            sample_1 = self.run_process("upload-cxb", inputs)

            inputs = {
                "src": "cuffquant_2.cxb",
                "source": "UCSC",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            sample_2 = self.run_process("upload-cxb", inputs)

            inputs = {
                "src": "3-cuffquant.cxb",
                "source": "UCSC",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            sample_3 = self.run_process("upload-cxb", inputs)

            inputs = {
                "src": "4-cuffquant.cxb",
                "source": "UCSC",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            sample_4 = self.run_process("upload-cxb", inputs)

            inputs = {
                "src": "5-cuffquant.cxb",
                "source": "UCSC",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            sample_5 = self.run_process("upload-cxb", inputs)

            inputs = {
                "src": "6-cuffquant.cxb",
                "source": "UCSC",
                "species": "Homo sapiens",
                "build": "hg19",
            }
            sample_6 = self.run_process("upload-cxb", inputs)

            RelationPartition.objects.create(
                relation=replicate_group, entity=sample_1.entity, label="1"
            )
            RelationPartition.objects.create(
                relation=replicate_group, entity=sample_2.entity, label="1"
            )
            RelationPartition.objects.create(
                relation=replicate_group, entity=sample_3.entity, label="2"
            )
            RelationPartition.objects.create(
                relation=replicate_group, entity=sample_4.entity, label="2"
            )
            RelationPartition.objects.create(
                relation=replicate_group, entity=sample_5.entity, label="2"
            )
            RelationPartition.objects.create(
                relation=replicate_group, entity=sample_6.entity, label="3"
            )

            annotation = self.prepare_annotation(
                fn="hg19_chr20_small.gtf.gz",
                source="UCSC",
                species="Homo sapiens",
                build="hg19",
            )

            self.assertEqual(
                replicate_groups(
                    [
                        {"__id": sample_1.id},
                        {"__id": sample_2.id},
                        {"__id": sample_3.id},
                        {"__id": sample_4.id},
                        {"__id": sample_5.id},
                        {"__id": sample_6.id},
                    ]
                ),
                [1, 1, 2, 2, 2, 3],
            )

        inputs = {
            "cuffquant": [
                sample_1.pk,
                sample_2.pk,
                sample_3.pk,
                sample_4.pk,
                sample_5.pk,
                sample_6.pk,
            ],
            "annotation": annotation.id,
        }
        cuffnorm = self.run_process("cuffnorm", inputs)
        self.assertFile(cuffnorm, "fpkm_means", "cuffnorm_all_fpkm_means.txt")
        self.assertFile(cuffnorm, "genes_fpkm", "cuffnorm_genes.fpkm_table")
        self.assertFileExists(cuffnorm, "raw_scatter")
        self.assertFields(cuffnorm, "source", "UCSC")
        self.assertFields(cuffnorm, "species", "Homo sapiens")
        self.assertFields(cuffnorm, "build", "hg19")

        exp = Data.objects.last()
        self.assertFile(exp, "exp", "cuffnorm_expression.tab.gz", compression="gzip")
        self.assertFile(
            exp, "exp_set", "cuffnorm_out_exp_set.txt.gz", compression="gzip"
        )
        self.assertJSON(exp, exp.output["exp_set_json"], "", "cuffnorm_exp_set.json.gz")

    @tag_process("mappability-bcm")
    def test_mappability(self):
        with self.preparation_stage():
            genome = self.prepare_ref_seq("genome.fasta.gz")
            bowtie_index = self.run_process("bowtie-index", {"ref_seq": genome.id})
            annotation = self.prepare_annotation_gff()

        mappability = self.run_process(
            "mappability-bcm",
            {
                "genome": bowtie_index.id,
                "gff": annotation.id,
                "length": 50,
            },
        )

        self.assertFileExists(mappability, "mappability")

    @tag_process("expression-dicty")
    def test_expression_dicty(self):
        with self.preparation_stage():
            ref_seq = self.prepare_ref_seq(
                fn="genome.fasta.gz",
                species="Dictyostelium discoideum",
                build="dd-05-2009",
            )
            hisat2_index = self.run_process("hisat2-index", {"ref_seq": ref_seq.id})
            reads = self.prepare_reads()
            annotation = self.prepare_annotation_gff()

            aligned_reads = self.run_process(
                "alignment-hisat2", {"genome": hisat2_index.pk, "reads": reads.pk}
            )

            mappa = self.run_process(
                "upload-mappability", {"src": "purpureum_mappability_50.tab.gz"}
            )

        inputs = {
            "alignment": aligned_reads.pk,
            "gff": annotation.pk,
            "mappable": mappa.pk,
        }
        expression = self.run_process("expression-dicty", inputs)
        self.assertFile(
            expression, "rpkm", "expression_bcm_rpkm.tab.gz", compression="gzip"
        )
        self.assertFields(expression, "source", "DICTYBASE")
        self.assertFields(expression, "species", "Dictyostelium discoideum")
        self.assertFields(expression, "build", "dd-05-2009")
        self.assertFields(expression, "feature_type", "gene")

    @with_resolwe_host
    @tag_process("mergeexpressions")
    def test_mergeexpression(self):
        with self.preparation_stage():
            expression_1 = self.prepare_expression(
                f_rc="exp_1_rc.tab.gz", f_exp="exp_1_tpm.tab.gz", f_type="TPM"
            )
            expression_2 = self.prepare_expression(
                f_rc="exp_2_rc.tab.gz", f_exp="exp_2_tpm.tab.gz", f_type="TPM"
            )
            expression_3 = self.prepare_expression(
                f_rc="exp_2_rc.tab.gz", f_exp="exp_2_tpm.tab.gz", f_type="RC"
            )

        inputs = {
            "exps": [expression_1.pk, expression_2.pk],
            "genes": ["DPU_G0067096", "DPU_G0067098", "DPU_G0067102"],
        }

        mergeexpression_1 = self.run_process("mergeexpressions", inputs)
        self.assertFile(mergeexpression_1, "expset", "merged_expset_subset.tab")

        inputs = {"exps": [expression_1.pk, expression_2.pk], "genes": []}

        mergeexpression_2 = self.run_process("mergeexpressions", inputs)
        self.assertFile(mergeexpression_2, "expset", "merged_expset_all.tab")

        inputs = {
            "exps": [expression_1.pk, expression_2.pk, expression_3.pk],
            "genes": ["DPU_G0067096", "DPU_G0067098", "DPU_G0067102"],
        }
        self.run_process("mergeexpressions", inputs, Data.STATUS_ERROR)

    @with_resolwe_host
    @tag_process("feature_counts")
    def test_feature_counts(self):
        base = Path("test_featurecounts")
        inputs = base / "inputs"
        outputs = base / "outputs"
        input_folder = Path("test_star") / "input"
        with self.preparation_stage():
            annotation_gtf = self.run_process(
                "upload-gtf",
                {
                    "src": inputs / "feature_counts hs.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens90",
                },
            )
            annotation_gff3 = self.prepare_annotation_gff()

            annotation = self.prepare_annotation(
                fn=input_folder / "hs annotation.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="GRCh38_ens90",
            )

            bam_single = self.run_process(
                "upload-bam",
                {
                    "src": inputs / "reads.bam",
                    "species": "Dictyostelium discoideum",
                    "build": "dd-05-2009",
                },
            )

            single_lanes = self.run_process(
                "upload-bam",
                {
                    "src": inputs / "single_lanes.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens90",
                },
            )

            bam_paired = self.run_process(
                "upload-bam",
                {
                    "src": inputs / "feature_counts hs_paired.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens90",
                },
            )

            paired_lanes = self.run_process(
                "upload-bam",
                {
                    "src": inputs / "paired_lanes.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens90",
                },
            )

            bam_ucsc = self.run_process(
                "upload-bam",
                {
                    "src": inputs / "cuffquant_mapping.bam",
                    "species": "Homo sapiens",
                    "build": "hg19",
                },
            )

            annotation_ucsc = self.prepare_annotation(
                fn=inputs / "hg19_chr20_small_modified.gtf.gz",
                source="UCSC",
                species="Homo sapiens",
                build="hg19",
            )
        # test using BAM file containing paired-end reads and a GTF input file
        expression = self.run_process(
            "feature_counts",
            {
                "aligned_reads": bam_paired.id,
                "annotation": annotation_gtf.id,
            },
        )
        self.assertFile(
            obj=expression,
            field_path="rc",
            fn=outputs / "feature_counts_out_rc.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            obj=expression,
            field_path="cpm",
            fn=outputs / "feature_counts_out_cpm.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            obj=expression,
            field_path="exp",
            fn=outputs / "feature_counts_out_tpm.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            obj=expression,
            field_path="exp_set",
            fn=outputs / "feature_counts_out_exp_set.txt.gz",
            compression="gzip",
        )
        self.assertJSON(
            expression,
            expression.output["exp_set_json"],
            "",
            outputs / "feature_counts_exp_set.json.gz",
        )
        self.assertFields(expression, "species", "Homo sapiens")
        self.assertFields(expression, "build", "GRCh38_ens90")
        self.assertFields(expression, "feature_type", "gene")

        # test using BAM file containing single-end reads and a GFF input file
        expression = self.run_process(
            "feature_counts",
            {
                "aligned_reads": bam_single.id,
                "annotation": annotation_gff3.id,
                "id_attribute": "Parent",
            },
        )
        self.assertFile(
            obj=expression,
            field_path="rc",
            fn=outputs / "feature_counts_out_gff3_rc.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            obj=expression,
            field_path="exp",
            fn=outputs / "feature_counts_out_gff3_tpm.tab.gz",
            compression="gzip",
        )
        self.assertFields(expression, "feature_type", "gene")

        expression_lanes = self.run_process(
            "feature_counts",
            {
                "aligned_reads": single_lanes.id,
                "annotation": annotation.id,
            },
        )
        self.assertFile(
            obj=expression_lanes,
            field_path="rc",
            fn=outputs / "feature_counts_single_lanes_rc.tab.gz",
            compression="gzip",
        )
        # test using UCSC-derived annotation
        expression = self.run_process(
            "feature_counts",
            {
                "aligned_reads": bam_ucsc.id,
                "annotation": annotation_ucsc.id,
            },
        )
        self.assertFile(
            obj=expression,
            field_path="rc",
            fn=outputs / "feature_counts_out_ucsc_rc.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            obj=expression,
            field_path="exp",
            fn=outputs / "feature_counts_out_ucsc_tpm.tab.gz",
            compression="gzip",
        )

        expression_lanes = self.run_process(
            "feature_counts",
            {
                "aligned_reads": paired_lanes.id,
                "annotation": annotation.id,
            },
        )
        self.assertFile(
            obj=expression_lanes,
            field_path="exp",
            fn=outputs / "feature_counts_paired_lanes_tpm.tab.gz",
            compression="gzip",
        )

    @tag_process("salmon-index")
    def test_salmon_index(self):
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
            "gencode": False,
            "keep_duplicates": True,
            "source": "ENSEMBL",
            "species": "Homo sapiens",
            "build": "ens_90",
        }
        salmon_index = self.run_process("salmon-index", inputs)

        del salmon_index.output["index"]["total_size"]  # Non-deterministic output.
        self.assertFields(salmon_index, "index", {"dir": "salmon_index"})
        self.assertFields(salmon_index, "source", "ENSEMBL")
        self.assertFields(salmon_index, "species", "Homo sapiens")
        self.assertFields(salmon_index, "build", "ens_90")

    @with_resolwe_host
    @tag_process("salmon-quant")
    def test_salmon_quant(self):
        inputs = Path("salmon_quant", "input")
        outputs = Path("salmon_quant", "output")

        with self.preparation_stage():
            reads = self.prepare_reads(
                [os.path.join("salmon_quant", "input", "hs sim_reads_single.fastq.gz")]
            )
            annotation = self.prepare_annotation(
                os.path.join("salmon_quant", "input", "hs annotation.gtf.gz"),
                source="ENSEMBL",
                species="Homo sapiens",
                build="ens_92",
            )
            transcripts = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": os.path.join("salmon_quant", "input", "hs cdna.fasta.gz"),
                    "species": "Homo sapiens",
                    "build": "ens_92",
                },
            )
            salmon_index = self.run_process(
                "salmon-index",
                {
                    "nucl": transcripts.id,
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens_92",
                },
            )

        inputs = {
            "reads": reads.id,
            "salmon_index": salmon_index.id,
            "annotation": annotation.id,
            "options": {
                "min_assigned_frag": 5,
                "gc_bias": True,
                "seq_bias": True,
                "incompat_prior": 0.05,
                "min_score_fraction": 0.7,
                "consensus_slack": 0.25,
                "no_length_correction": False,
                "discard_orphans_quasi": True,
            },
        }
        salmon_quant = self.run_process("salmon-quant", inputs)
        self.assertFile(
            obj=salmon_quant,
            field_path="exp_set",
            fn=outputs / "salmon_quant_tpm.tab.gz",
            compression="gzip",
        )
        self.assertFile(
            obj=salmon_quant,
            field_path="transcripts",
            fn=outputs / "salmon_transcripts_tpm.tab.gz",
            compression="gzip",
        )

        self.assertFile(
            obj=salmon_quant,
            field_path="rc",
            fn=outputs / "salmon_counts.txt.gz",
            compression="gzip",
        )

        inputs = {
            "reads": reads.id,
            "salmon_index": salmon_index.id,
            "annotation": annotation.id,
            "options": {
                "min_assigned_frag": 5,
                "gc_bias": True,
                "seq_bias": True,
                "incompat_prior": 0.05,
                "min_score_fraction": 0.7,
                "consensus_slack": 0.25,
                "no_length_correction": False,
                "discard_orphans_quasi": True,
                "num_bootstraps": 5,
            },
        }
        salmon_quant = self.run_process("salmon-quant", inputs)
        self.assertFile(
            obj=salmon_quant,
            field_path="transcripts",
            fn=outputs / "salmon_transcripts_tpm.tab.gz",
            compression="gzip",
        )

        self.assertFile(
            obj=salmon_quant,
            field_path="variance",
            fn=outputs / "variance_salmon.txt.gz",
            compression="gzip",
        )

    @with_resolwe_host
    @tag_process("feature_counts")
    def test_featurecounts_strandedness(self):
        base = Path("test_featurecounts")
        inputs = base / "inputs"
        outputs = base / "outputs"
        with self.preparation_stage():
            cds = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": str(inputs / "salmon_cds.fa.gz"),
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            salmon_index = self.run_process(
                "salmon-index",
                {
                    "nucl": cds.id,
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            annotation = self.run_process(
                "upload-gtf",
                {
                    "src": str(inputs / "annotation_rsem.gtf.gz"),
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            aligned_reads_paired = self.run_process(
                "upload-bam",
                {
                    "src": "feature counts_detect_strandedness.bam",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

            aligned_reads_single = self.run_process(
                "upload-bam",
                {
                    "src": "feature counts_detect_strandedness_single.bam",
                    "species": "Homo sapiens",
                    "build": "ens_90",
                },
            )

        expression_paired = self.run_process(
            "feature_counts",
            {
                "aligned_reads": aligned_reads_paired.id,
                "assay_type": "auto",
                "cdna_index": salmon_index.id,
                "annotation": annotation.id,
            },
        )
        self.assertFile(
            obj=expression_paired,
            field_path="exp",
            fn=outputs / "auto_detect_strand_tpm.tab.gz",
            compression="gzip",
        )

        expression_single = self.run_process(
            "feature_counts",
            {
                "aligned_reads": aligned_reads_single.id,
                "assay_type": "auto",
                "cdna_index": salmon_index.id,
                "annotation": annotation.id,
            },
        )
        self.assertFile(
            obj=expression_single,
            field_path="exp",
            fn=outputs / "auto_detect_strand_tpm.tab.gz",
            compression="gzip",
        )

    @with_resolwe_host
    @tag_process("star-quantification")
    def test_star_quantification(self):
        inputs = Path("test_star") / "input"
        outputs = Path("test_star") / "output"
        with self.preparation_stage():
            paired_lanes = self.prepare_paired_reads(
                mate1=[
                    inputs / "hs_paired_R1 workflow_bbduk_star_htseq.fastq.gz",
                    "hs sim_reads1.fastq.gz",
                ],
                mate2=[
                    inputs / "hs_paired_R2 workflow_bbduk_star_htseq.fastq.gz",
                    "hs sim_reads2.fastq.gz",
                ],
            )
            single_lanes = self.prepare_reads(
                [
                    inputs / "hs_single bbduk_star_htseq_reads_single.fastq.gz",
                    "hs sim_reads_single.fastq.gz",
                ]
            )
            annotation = self.prepare_annotation(
                fn=inputs / "hs annotation.gtf.gz",
                source="ENSEMBL",
                species="Homo sapiens",
                build="GRCh38_ens90",
            )
            inputs = {
                "src": inputs / "hs genome.fasta.gz",
                "species": "Homo sapiens",
                "build": "GRCh38_ens90",
            }
            star_index_fasta = self.run_process("upload-fasta-nucl", inputs)

            star_index = self.run_process(
                "alignment-star-index",
                {"annotation": annotation.id, "ref_seq": star_index_fasta.id},
            )

            salmon_index = self.run_process(
                "salmon-index",
                {
                    "nucl": star_index_fasta.id,
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens90",
                },
            )

            inputs = {
                "genome": star_index.id,
                "reads": paired_lanes.id,
                "gene_counts": True,
            }
            aligned_reads = self.run_process("alignment-star", inputs)

            inputs = {
                "genome": star_index.id,
                "reads": single_lanes.id,
                "gene_counts": True,
            }
            aligned_reads_single = self.run_process("alignment-star", inputs)

            inputs = {
                "genome": star_index.id,
                "reads": single_lanes.id,
            }
            aligned_reads_error = self.run_process("alignment-star", inputs)

        # test using BAM file containing paired-end reads and a GTF input file
        expression = self.run_process(
            "star-quantification",
            {
                "aligned_reads": aligned_reads.id,
                "annotation": annotation.id,
            },
        )
        self.assertFields(expression, "species", "Homo sapiens")
        self.assertFields(expression, "build", "GRCh38_ens90")
        self.assertFile(
            obj=expression,
            field_path="exp",
            fn=outputs / "star_tpm.tab.gz",
            compression="gzip",
        )
        expression = self.run_process(
            "star-quantification",
            {
                "aligned_reads": aligned_reads_single.id,
                "annotation": annotation.id,
            },
        )
        self.assertFields(expression, "species", "Homo sapiens")
        self.assertFields(expression, "build", "GRCh38_ens90")
        self.assertFile(
            obj=expression,
            field_path="exp_set",
            fn=outputs / "star_single_expressions.txt.gz",
            compression="gzip",
        )

        expression = self.run_process(
            "star-quantification",
            {
                "aligned_reads": aligned_reads_single.id,
                "annotation": annotation.id,
                "assay_type": "auto",
                "cdna_index": salmon_index.id,
            },
        )
        self.assertFields(expression, "species", "Homo sapiens")
        self.assertFields(expression, "build", "GRCh38_ens90")
        self.assertFile(
            expression,
            "counts_summary",
            outputs / "exp_summary.txt",
        )
        self.assertFields(expression, "feature_type", "gene")
        expression = self.run_process(
            "star-quantification",
            {
                "aligned_reads": aligned_reads_error.id,
                "annotation": annotation.id,
                "assay_type": "auto",
                "cdna_index": salmon_index.id,
            },
            Data.STATUS_ERROR,
        )
        self.assertEqual(
            expression.process_error,
            ["Aligned reads should contain gene count information, but do not."],
        )
