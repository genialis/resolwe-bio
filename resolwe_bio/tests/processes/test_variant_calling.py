from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.filter import filter_html, filter_vcf_variable
from resolwe_bio.utils.test import BioProcessTestCase, skipUnlessLargeFiles


class VariantCallingTestCase(BioProcessTestCase):
    @tag_process("vc-chemut")
    def test_variant_calling_chemut(self):
        input_folder = Path("chemut") / "input"
        output_folder = Path("chemut") / "output"
        with self.preparation_stage():
            inputs = {
                "src": input_folder / "chemut_genome.fasta.gz",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            ref_seq = self.run_process("upload-fasta-nucl", inputs)
            bwa_index = self.run_process("bwa-index", {"ref_seq": ref_seq.id})

            inputs = {
                "src1": [input_folder / "AX4_mate1.fq.gz"],
                "src2": [input_folder / "AX4_mate2.fq.gz"],
            }

            parental_reads = self.run_process("upload-fastq-paired", inputs)

            inputs = {
                "src1": [input_folder / "CM_mate1.fq.gz"],
                "src2": [input_folder / "CM_mate2.fq.gz"],
            }

            mut_reads = self.run_process("upload-fastq-paired", inputs)

            inputs = {"genome": bwa_index.id, "reads": parental_reads.id}
            align_parental = self.run_process("alignment-bwa-mem", inputs)

            inputs = {"genome": bwa_index.id, "reads": mut_reads.id}
            align_mut = self.run_process("alignment-bwa-mem", inputs)

        inputs = {
            "genome": ref_seq.id,
            "parental_strains": [align_parental.id],
            "mutant_strains": [align_mut.id],
            "reads_info": {
                "PL": "Illumina",
                "LB": "x",
                "CN": "def",
                "DT": "2014-08-05",
            },
        }

        variants = self.run_process("vc-chemut", inputs)
        self.assertFile(
            variants,
            "vcf",
            output_folder / "GATKvariants_raw.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )
        self.assertFields(variants, "build", "dd-05-2009")
        self.assertFields(variants, "species", "Dictyostelium discoideum")

    @tag_process("filtering-chemut")
    def test_filtering_chemut(self):
        input_folder = Path("chemut") / "input"
        output_folder = Path("chemut") / "output"
        with self.preparation_stage():
            vcf_input = {
                "src": input_folder / "variant_calling_filtering.vcf.gz",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            variants = self.run_process("upload-variants-vcf", vcf_input)
            inputs = {
                "src": input_folder / "dd_chr1.fasta.gz",
                "species": "Dictyostelium discoideum",
                "build": "dd-05-2009",
            }
            ref_seq = self.run_process("upload-fasta-nucl", inputs)

        inputs = {
            "variants": variants.pk,
            "analysis_type": "snv",
            "parental_strain": "parental",
            "mutant_strain": "mut",
            "genome": ref_seq.id,
            "read_depth": 5,
        }

        filtered_variants = self.run_process("filtering-chemut", inputs)
        self.assertFile(
            filtered_variants,
            "vcf",
            output_folder / "snv_filtered.vcf.gz",
            compression="gzip",
        )
        self.assertFile(
            filtered_variants,
            "variants_filtered",
            output_folder / "variant_filtered.txt",
        )
        self.assertFields(filtered_variants, "build", "dd-05-2009")
        self.assertFields(filtered_variants, "species", "Dictyostelium discoideum")

        inputs = {
            "variants": variants.pk,
            "analysis_type": "indel",
            "parental_strain": "parental",
            "mutant_strain": "mut",
            "genome": ref_seq.id,
            "read_depth": 5,
        }

        filtered_variants = self.run_process("filtering-chemut", inputs)
        self.assertFile(
            filtered_variants,
            "vcf",
            output_folder / "indel_filtered.vcf.gz",
            compression="gzip",
        )
        self.assertFields(filtered_variants, "build", "dd-05-2009")
        self.assertFields(filtered_variants, "species", "Dictyostelium discoideum")

    @skipUnlessLargeFiles("56GSID_10k.realigned.bqsrCal.bam")
    @tag_process("vc-gatk4-hc")
    def test_gatk_haplotypecaller(self):
        with self.preparation_stage():
            input_folder = Path("haplotypecaller") / "input"
            output_folder = Path("haplotypecaller") / "output"
            alignment = self.run_process(
                "upload-bam",
                {
                    "src": Path("large") / "56GSID_10k.realigned.bqsrCal.bam",
                    "species": "Homo sapiens",
                    "build": "b37",
                },
            )

            inputs = {
                "src": "hs_b37_chr2_small.fasta.gz",
                "species": "Homo sapiens",
                "build": "b37",
            }
            genome = self.run_process("upload-fasta-nucl", inputs)

            bed_file = self.run_process(
                "upload-bed",
                {
                    "src": input_folder / "56G_masterfile_test_merged_targets_5col.bed",
                    "species": "Homo sapiens",
                    "build": "b37",
                },
            )

            dbsnp_input = {
                "src": "dbsnp_138.b37.chr2_small.vcf.gz",
                "species": "Homo sapiens",
                "build": "b37",
            }
            dbsnp = self.run_process("upload-variants-vcf", dbsnp_input)

        gatk4 = self.run_process(
            "vc-gatk4-hc",
            {
                "alignment": alignment.id,
                "intervals_bed": bed_file.id,
                "genome": genome.id,
                "dbsnp": dbsnp.id,
            },
        )

        self.assertFile(
            gatk4,
            "vcf",
            output_folder / "56GSID_10k.gatkHC4.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

        # Test specifically for HaplotypeCaller with RNA-seq data
        gatk_rnaseq = self.run_process(
            "vc-gatk4-hc",
            {
                "alignment": alignment.id,
                "genome": genome.id,
                "dbsnp": dbsnp.id,
                "mbq": 10,
                "advanced": {
                    "soft_clipped": True,
                    "java_gc_threads": 3,
                    "max_heap_size": 7,
                },
            },
        )
        self.assertFields(gatk_rnaseq, "build", "b37")
        self.assertFields(gatk_rnaseq, "species", "Homo sapiens")
        self.assertFile(
            gatk_rnaseq,
            "vcf",
            output_folder / "56GSID_10k.rna-seq.gatkHC.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

    @tag_process("snpeff", "snpeff-single")
    def test_snpeff(self):
        with self.preparation_stage():
            input_folder = Path("snpeff") / "input"
            output_folder = Path("snpeff") / "output"
            multi = Path("wgs") / "input"
            variants_rna = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "filtered_snpeff.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38_ens100",
                },
            )
            dbsnp_rna = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "dbsnp.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh37",
                },
            )
            dbsnp_rna_38 = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "dbsnp.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            genes = self.run_process(
                "upload-geneset",
                {
                    "src": input_folder / "set_of_genes.txt",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                },
            )
            vcf_multi = self.run_process(
                "upload-variants-vcf",
                {
                    "src": multi / "vcf_1.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

        snpeff = self.run_process(
            "snpeff-single",
            {
                "variants": variants_rna.id,
                "database": "GRCh38.109",
                "dbsnp": dbsnp_rna.id,
            },
            Data.STATUS_ERROR,
        )
        self.assertEqual(
            snpeff.process_error[0],
            "Genome build for the DBSNP file and used database "
            "should be the same. DBSNP file is based on "
            "GRCh37, while snpEff database is based on "
            "GRCh38.",
        )

        snpeff = self.run_process(
            "snpeff-single",
            {
                "variants": vcf_multi.id,
                "database": "GRCh38.109",
            },
            Data.STATUS_ERROR,
        )
        self.assertEqual(
            snpeff.process_error[0],
            "The input VCF should contain data for a single sample. "
            "The input contains data for 2 sample(s).",
        )

        snpeff = self.run_process(
            "snpeff",
            {
                "variants": variants_rna.id,
                "database": "GRCh38.109",
            },
            Data.STATUS_ERROR,
        )
        self.assertEqual(
            snpeff.process_error[0],
            "The input VCF file should contain data for multiple samples. "
            "The input contains data for 1 sample(s).",
        )

        snpeff = self.run_process(
            "snpeff-single",
            {
                "variants": variants_rna.id,
                "database": "GRCh38.109",
                "dbsnp": dbsnp_rna_38.id,
            },
        )
        self.assertFields(snpeff, "build", "GRCh38_ens100")

        snpeff_filtering = self.run_process(
            "snpeff-single",
            {
                "variants": variants_rna.id,
                "database": "GRCh38.109",
                "filtering_options": "( REF = 'A' )",
                "extract_fields": [
                    "CHROM",
                    "POS",
                    "REF",
                    "ALT",
                    "ANN[*].GENE",
                    "ANN[*].HGVS_P",
                ],
                "advanced": {"one_per_line": True},
            },
        )
        self.assertFile(
            snpeff_filtering,
            "vcf_extracted",
            output_folder / "extracted_variants.vcf.gz",
            compression="gzip",
            file_filter=filter_vcf_variable,
        )
        self.assertFile(
            snpeff_filtering,
            "vcf",
            output_folder / "filtered_variants.vcf.gz",
            compression="gzip",
            file_filter=filter_vcf_variable,
        )

        snpeff_filtering = self.run_process(
            "snpeff-single",
            {
                "variants": variants_rna.id,
                "database": "GRCh38.109",
                "filtering_options": "ANN[*].EFFECT has 'missense_variant'",
                "sets": [genes.id],
                "extract_fields": [
                    "CHROM",
                    "POS",
                    "ID",
                    "REF",
                    "ALT",
                    "QUAL",
                    "ANN[*].GENE",
                    "ANN[*].HGVS_P",
                ],
            },
        )
        self.assertFile(
            snpeff_filtering,
            "vcf_extracted",
            output_folder / "variants_set.vcf.gz",
            compression="gzip",
            file_filter=filter_vcf_variable,
        )

        snpeff = self.run_process(
            "snpeff",
            {
                "variants": vcf_multi.id,
                "database": "GRCh38.109",
            },
        )
        self.assertFile(
            snpeff,
            "vcf",
            output_folder / "multi_variants.vcf.gz",
            compression="gzip",
            file_filter=filter_vcf_variable,
        )

    @tag_process("gatk-refine-variants")
    def test_gatk_refinement(self):
        input_folder = Path("variant_refinement") / "input"
        output_folder = Path("variant_refinement") / "output"
        with self.preparation_stage():
            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": input_folder / "chrX_1_28000.fa.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )
            variants_main = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "snp.recalibrated_chrX_28000.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )
            vcf_pop = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "chrX_28000.renamed.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

        gatk_cgp = self.run_process(
            "gatk-refine-variants",
            {
                "vcf": variants_main.id,
                "ref_seq": ref_seq.id,
                "vcf_pop": vcf_pop.id,
            },
        )
        self.assertFile(
            gatk_cgp,
            "vcf",
            output_folder / "variants.refined.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

    @tag_process("ensembl-vep")
    def test_ensembl_vep(self):
        input_folder = Path("ensembl-vep") / "input"
        output_folder = Path("ensembl-vep") / "output"
        with self.preparation_stage():
            vcf = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "snp.recalibrated_chrX_28000.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            cache = self.run_process(
                "upload-vep-cache",
                {
                    "cache_file": input_folder / "cache_homo_sapiens_X.tar.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                    "release": "103",
                },
            )
            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": input_folder / "chrX_1_28000.fa.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )
        vep = self.run_process(
            "ensembl-vep",
            {
                "vcf": vcf.id,
                "cache": cache.id,
                "ref_seq": ref_seq.id,
            },
        )
        self.assertFile(
            vep,
            "vcf",
            output_folder / "snp_annotated.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )
        self.assertFile(vep, "tbi", output_folder / "snp_annotated.vcf.gz.tbi")
        self.assertFile(
            vep,
            "summary",
            output_folder / "snp_annotated.vcf_summary.html",
            file_filter=filter_html,
        )
        self.assertEqual(
            vep.process_warning,
            [
                "The current version of Ensembl-VEP is 104. "
                "It is recommended that cache version is also 104."
            ],
        )

    @tag_process("variants-to-table")
    def test_variants_to_table(self):
        input_folder = Path("variants_to_table") / "input"
        output_folder = Path("variants_to_table") / "output"
        with self.preparation_stage():
            vcf = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "refined_variants.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
        variants = self.run_process(
            "variants-to-table",
            {
                "vcf": vcf.id,
                "vcf_fields": ["CHROM", "POS", "ALT"],
                "advanced_options": {
                    "gf_fields": ["GT", "GQ"],
                    "split_alleles": False,
                },
            },
        )
        self.assertFile(variants, "tsv", output_folder / "variants_to_table.tsv")

    @tag_process("gatk-variant-filtration", "gatk-variant-filtration-single")
    def test_variant_filtration(self):
        input_folder = Path("variant_filtration") / "input"
        output_folder = Path("variant_filtration") / "output"
        input_dbsnp = Path("rnaseq_variantcalling") / "input"
        with self.preparation_stage():
            ref_seq = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": input_folder / "chr1_19000.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )
            vcf = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "chr1_19000.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            vcf_multi_sample = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "chr1_19000_multi_sample.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "custom_build",
                },
            )

            dbsnp = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_dbsnp / "dbsnp-hg38.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

        filtering_single_sample = self.run_process(
            "gatk-variant-filtration-single",
            {
                "vcf": vcf.id,
                "ref_seq": ref_seq.id,
                "filter_expressions": ["FS > 30.0", "QD < 2.0", "DP > 20"],
                "filter_name": ["FS", "QD", "DP"],
                "advanced": {
                    "window": 35,
                    "java_gc_threads": 3,
                    "max_heap_size": 8,
                },
            },
        )
        self.assertFile(
            filtering_single_sample,
            "vcf",
            output_folder / "filtered_variants.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

        filtering_mask = self.run_process(
            "gatk-variant-filtration-single",
            {
                "vcf": vcf.id,
                "ref_seq": ref_seq.id,
                "filter_expressions": ["FS > 30.0", "QD < 2.0", "DP > 20"],
                "filter_name": ["FS", "QD", "DP"],
                "mask": dbsnp.id,
                "mask_name": "DB",
                "advanced": {
                    "window": 35,
                    "java_gc_threads": 3,
                    "max_heap_size": 8,
                },
            },
        )
        self.assertFile(
            filtering_mask,
            "vcf",
            output_folder / "filtered_variants_mask.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

        filtering_single_sample = self.run_process(
            "gatk-variant-filtration-single",
            {
                "vcf": vcf.id,
                "ref_seq": ref_seq.id,
                "genotype_filter_expressions": [
                    "DP < 8.0",
                    "vc.getGenotype('20').getAD().1 < 3.0",
                ],
                "genotype_filter_name": ["DP", "AD"],
            },
        )
        self.assertFile(
            filtering_single_sample,
            "vcf",
            output_folder / "genotype_filtered_variants.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

        filtering_error = {
            "vcf": vcf.id,
            "ref_seq": ref_seq.id,
            "filter_expressions": ["QD < 2.0", "DP > 20"],
            "filter_name": ["FS", "QD", "DP"],
        }

        filtering = self.run_process(
            "gatk-variant-filtration-single", filtering_error, Data.STATUS_ERROR
        )
        error_msg = [
            ("The number of filter expressions and filter names is not the same.")
        ]
        self.assertEqual(filtering.process_error, error_msg)

        filtering_multi_sample = self.run_process(
            "gatk-variant-filtration",
            {
                "vcf": vcf_multi_sample.id,
                "ref_seq": ref_seq.id,
                "filter_expressions": ["FS > 30.0", "QD < 2.0", "DP > 20"],
                "filter_name": ["FS", "QD", "DP"],
                "advanced": {
                    "window": 35,
                    "java_gc_threads": 3,
                    "max_heap_size": 8,
                },
            },
        )
        self.assertFile(
            filtering_multi_sample,
            "vcf",
            output_folder / "filtered_variants_multi_sample.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

        filtering_multi_sample = self.run_process(
            "gatk-variant-filtration",
            {
                "vcf": vcf_multi_sample.id,
                "ref_seq": ref_seq.id,
                "genotype_filter_expressions": ["AD.1 < 5"],
                "genotype_filter_name": ["AD"],
            },
        )
        self.assertFile(
            filtering_multi_sample,
            "vcf",
            output_folder / "genotype_filtered_variants_multi.vcf.gz",
            file_filter=filter_vcf_variable,
            compression="gzip",
        )

        wrong_process_inputs = {
            "vcf": vcf_multi_sample.id,
            "ref_seq": ref_seq.id,
            "filter_expressions": ["FS > 30.0", "QD < 2.0", "DP > 20"],
            "filter_name": ["FS", "QD", "DP"],
            "advanced": {
                "window": 35,
                "java_gc_threads": 3,
                "max_heap_size": 8,
            },
        }

        wrong_process = self.run_process(
            "gatk-variant-filtration-single", wrong_process_inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "The input VCF should contain data for a single sample. "
            "The input contains data for 2 sample(s)."
        ]
        self.assertEqual(wrong_process.process_error, error_msg)
