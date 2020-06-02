import os

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class UploadProcessorTestCase(KBBioProcessTestCase):
    @tag_process("upload-bam", "upload-bam-indexed")
    def test_bam_upload(self):
        inputs = {
            "src": "alignment_name_sorted.bam",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        upload_bam = self.run_process("upload-bam", inputs)
        self.assertFile(upload_bam, "bam", "alignment_position_sorted.bam")
        self.assertFile(upload_bam, "bai", "alignment_bam_upload_index.bai")
        self.assertFile(upload_bam, "stats", "alignment_bam_upload_stats.txt")
        self.assertFile(upload_bam, "bigwig", "alignment_bam_upload_bigwig.bw")
        self.assertFields(upload_bam, "species", "Homo sapiens")
        self.assertFields(upload_bam, "build", "hg19")

        inputs = {
            "src": "alignment_position_sorted.bam",
            "src2": "alignment_bam_upload_index.bam.bai",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        upload_bam = self.run_process("upload-bam-indexed", inputs, Data.STATUS_ERROR)
        self.assertEqual(
            upload_bam.process_error[0],
            "BAI should have the same name as BAM with .bai extension",
        )

        inputs = {
            "src": "alignment_position_sorted.bam",
            "src2": "alignment_position_sorted.bam.bai",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        upload_bam = self.run_process("upload-bam-indexed", inputs)
        self.assertFile(upload_bam, "bam", "alignment_position_sorted.bam")
        self.assertFile(upload_bam, "bai", "alignment_position_sorted.bam.bai")
        self.assertFile(upload_bam, "stats", "alignment_bam_upload_stats.txt")
        self.assertFile(upload_bam, "bigwig", "alignment_bam_upload_bigwig.bw")
        self.assertFields(upload_bam, "species", "Homo sapiens")
        self.assertFields(upload_bam, "build", "hg19")

    @with_resolwe_host
    @tag_process("upload-expression")
    def test_upload_expression(self):
        inputs = {
            "exp_type": "TPM",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        self.run_process("upload-expression", inputs, Data.STATUS_ERROR)

        inputs = {
            "exp": "exp_1_tpm.tab.gz",
            "rc": "exp_1_rc.tab.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        self.run_process("upload-expression", inputs, Data.STATUS_ERROR)

        inputs = {
            "rc": "exp_1_rc.tab.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_3 = self.run_process("upload-expression", inputs)
        self.assertFile(exp_3, "rc", "exp_1_rc.tab.gz")
        self.assertFile(exp_3, "exp", "exp_1_rc.tab.gz")
        self.assertJSON(exp_3, exp_3.output["exp_json"], "", "exp_1.json.gz")
        self.assertFields(exp_3, "species", "Homo sapiens")
        self.assertFields(exp_3, "build", "hg19")
        self.assertFields(exp_3, "feature_type", "gene")

        inputs = {
            "exp": "exp_1_tpm.tab.gz",
            "exp_type": "TPM",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_4 = self.run_process("upload-expression", inputs)
        self.assertFile(exp_4, "exp", "exp_1_tpm.tab.gz")

        inputs = {
            "rc": "exp_1_rc.tab.gz",
            "exp": "exp_1_tpm.tab.gz",
            "exp_type": "TPM",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_5 = self.run_process("upload-expression", inputs)
        self.assertFields(exp_5, "exp_type", "TPM")
        self.assertFile(exp_5, "exp", "exp_1_tpm.tab.gz")
        self.assertFile(exp_5, "rc", "exp_1_rc.tab.gz")
        self.assertJSON(exp_5, exp_5.output["exp_json"], "", "exp_1_norm.json.gz")
        self.assertJSON(
            exp_5, exp_5.output["exp_set_json"], "", "upload_exp_norm_set.json.gz"
        )

        inputs = {
            "rc": "exp_mac_line_ending.txt.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_6 = self.run_process("upload-expression", inputs)
        self.assertJSON(exp_6, exp_6.output["exp_json"], "", "exp.json.gz")

        inputs = {
            "rc": "exp_unix_line_ending.txt.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_7 = self.run_process("upload-expression", inputs)
        self.assertJSON(exp_7, exp_7.output["exp_json"], "", "exp.json.gz")

        inputs = {
            "rc": "exp_windows_line_ending.txt.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_8 = self.run_process("upload-expression", inputs)
        self.assertJSON(exp_8, exp_8.output["exp_json"], "", "exp.json.gz")

        # Check handling of numerical feature_ids in expression file
        inputs = {
            "rc": "mm_ncbi_exp.tab.gz",
            "exp_name": "Expression",
            "source": "NCBI",
            "species": "Mus musculus",
            "build": "hg19",
        }
        exp_9 = self.run_process("upload-expression", inputs)
        self.assertFile(exp_9, "exp_set", "mm_ncbi_exp_set.txt.gz", compression="gzip")

    @with_resolwe_host
    @tag_process("upload-cxb", "upload-expression-cuffnorm")
    def test_upload_cuffquant_expr(self):
        inputs = {
            "src": "cuffquant 1.cxb",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        cxb = self.run_process("upload-cxb", inputs)

        inputs = {"exp": "cuffquant_exp.tab", "cxb": cxb.id}
        exp = self.run_process("upload-expression-cuffnorm", inputs)
        self.assertFields(exp, "feature_type", "gene")

    @tag_process("upload-fastq-paired")
    def test_upload_paired_end_reads(self):
        inputs = {
            "src1": ["mate1_reordered.fastq.gz", "mate1_diff_num_reads.fastq.gz"],
            "src2": ["mate2_reordered.fastq.gz"],
        }
        wrong_mates = self.run_process("upload-fastq-paired", inputs, Data.STATUS_ERROR)
        error_msg = [
            "The number of mate-pair files in split-lane samples must match. 2 and 1 "
            "input files were given for the -fq and -fq2 inputs, respectively."
        ]
        self.assertEqual(wrong_mates.process_error, error_msg)

        inputs = {
            "src1": ["mate1_reordered.fastq.gz", "mate1_reordered.fastq.gz"],
            "src2": ["mate2_reordered.fastq.gz"],
        }
        wrong_mates2 = self.run_process(
            "upload-fastq-paired", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Non-unique input file names detected: ['mate1_reordered.fastq.gz']."
        ]
        self.assertEqual(wrong_mates2.process_error, error_msg)

        inputs = {
            "src1": ["mate1_diff_num_reads.fastq.gz"],
            "src2": ["mate2_diff_num_reads.fastq.gz"],
        }
        diff_numb_reads = self.run_process(
            "upload-fastq-paired", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Format error in mate-pairs mate1_diff_num_reads.fastq.gz and mate2_diff_num_reads.fastq.gz. "
            "Error in sequence file at unknown line: Reads are improperly paired. "
            "There are more reads in file 2 than in file 1."
        ]
        self.assertEqual(diff_numb_reads.process_error, error_msg)

        inputs = {
            "src1": ["mate1_reordered.fastq.gz"],
            "src2": ["mate2_reordered.fastq.gz"],
        }
        unordered_reads = self.run_process(
            "upload-fastq-paired", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Format error in mate-pairs mate1_reordered.fastq.gz and mate2_reordered.fastq.gz. "
            "Error in sequence file at unknown line: Reads are improperly paired. Read "
            "name 'read1/1 some text' in file 1 does not match 'read2/2' in file 2."
        ]
        self.assertEqual(unordered_reads.process_error, error_msg)

        inputs = {"src1": ["mate1.fastq.gz"], "src2": ["mate2.fastq.gz"]}
        self.run_process("upload-fastq-paired", inputs, Data.STATUS_ERROR)

        inputs = {"src1": ["rRNA forw.fastq.gz"], "src2": ["rRNA_rew.fastq.gz"]}
        reads = self.run_process("upload-fastq-paired", inputs)
        self.assertFiles(reads, "fastq", ["rRNA forw.fastq.gz"], compression="gzip")
        self.assertFiles(reads, "fastq2", ["rRNA_rew.fastq.gz"], compression="gzip")
        del reads.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        self.assertFields(
            reads,
            "fastqc_url",
            [
                {
                    "file": "fastqc/rRNA forw_fastqc/fastqc_report.html",
                    "refs": ["fastqc/rRNA forw_fastqc"],
                    "size": 343222,
                }
            ],
        )
        del reads.output["fastqc_url2"][0]["total_size"]  # Non-deterministic output.
        self.assertFields(
            reads,
            "fastqc_url2",
            [
                {
                    "file": "fastqc/rRNA_rew_fastqc/fastqc_report.html",
                    "refs": ["fastqc/rRNA_rew_fastqc"],
                    "size": 323297,
                }
            ],
        )

        merged_lanes = self.run_process(
            "upload-fastq-paired",
            {
                "src1": ["old_encoding.fastq.gz", "old_encoding1.fastq.gz"],
                "src2": ["old_encoding_R2.fastq.gz", "old_encoding1_R2.fastq.gz"],
                "merge_lanes": True,
            },
        )
        self.assertFiles(
            merged_lanes,
            "fastq",
            ["paired_end_merged_lanes_mate1.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            merged_lanes,
            "fastq2",
            ["paired_end_merged_lanes_mate2.fastq.gz"],
            compression="gzip",
        )

    @tag_process("upload-fastq-single")
    def test_upload_single_end_reads(self):
        empty_input = self.run_process(
            "upload-fastq-single", {"src": ["empty.fastq.gz"]}, Data.STATUS_ERROR
        )
        error_msg = ["Input file empty.fastq.gz contains no read sequences."]
        self.assertEqual(empty_input.process_error, error_msg)

        garbage_input = self.run_process(
            "upload-fastq-single", {"src": ["garbage.fastq.gz"]}, Data.STATUS_ERROR
        )
        error_msg = [
            "Error in file garbage.fastq.gz. Error in FASTQ file at line 1: Line expected "
            "to start with '@', but found 'S'"
        ]
        self.assertEqual(garbage_input.process_error, error_msg)

        garbage_input_2 = self.run_process(
            "upload-fastq-single", {"src": ["garbage2.fastq.gz"]}, Data.STATUS_ERROR
        )
        error_msg = [
            "Error in file garbage2.fastq.gz. Error in FASTQ file at line 3: Sequence descriptions "
            "don't match ('Some random content' != '+++'). The second sequence description must "
            "be either empty or equal to the first description."
        ]
        self.assertEqual(garbage_input_2.process_error, error_msg)

        missing_qual = self.run_process(
            "upload-fastq-single", {"src": ["missing_qual.fastq.gz"]}, Data.STATUS_ERROR
        )
        error_msg = [
            "Error in file missing_qual.fastq.gz. Error in FASTQ file at line 16: "
            "Premature end of file encountered. The incomplete final record was: "
            "'@read4/1\nGACAGGCCGTTTGAATGTTGACGGGATGTT\n+\n'"
        ]
        self.assertEqual(missing_qual.process_error, error_msg)

        inputs = {"src": ["mate1.fastq.gz"]}
        self.run_process("upload-fastq-single", inputs, Data.STATUS_ERROR)

        inputs = {"src": ["rRNA forw.fastq.gz", "rRNA_rew.fastq.gz"]}
        reads = self.run_process("upload-fastq-single", inputs)

        self.assertFiles(
            reads,
            "fastq",
            ["rRNA_forw_single.fastq.gz", "rRNA_rew.fastq.gz"],
            compression="gzip",
        )
        del reads.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        del reads.output["fastqc_url"][1]["total_size"]  # Non-deterministic output.
        self.assertFields(
            reads,
            "fastqc_url",
            [
                {
                    "file": "fastqc/rRNA forw_fastqc/fastqc_report.html",
                    "refs": ["fastqc/rRNA forw_fastqc"],
                    "size": 343222,
                },
                {
                    "file": "fastqc/rRNA_rew_fastqc/fastqc_report.html",
                    "refs": ["fastqc/rRNA_rew_fastqc"],
                    "size": 323297,
                },
            ],
        )

        merged_lanes = self.run_process(
            "upload-fastq-single",
            {"src": ["rRNA forw.fastq.gz", "rRNA_rew.fastq.gz"], "merge_lanes": True,},
        )
        self.assertFiles(
            merged_lanes,
            "fastq",
            ["merged_single_end_reads.fastq.gz"],
            compression="gzip",
        )

    @tag_process("upload-diffexp")
    def test_upload_de(self):
        inputs = {
            "src": "./diff_exp/input/deseq2_output.tab.gz",
            "source": "DICTYBASE",
            "gene_id": "gene_id",
            "logfc": "log2FoldChange",
            "fdr": "padj",
            "pvalue": "pvalue",
            "stat": "stat",
            "species": "Dictyostelium discoideum",
            "build": "dd-05-2009",
            "feature_type": "gene",
        }
        diff_exp = self.run_process("upload-diffexp", inputs)

        self.assertFile(diff_exp, "raw", "./diff_exp/input/deseq2_output.tab.gz")
        self.assertJSON(
            diff_exp,
            diff_exp.output["de_json"],
            "",
            "./diff_exp/output/deseq2_volcano_plot.json.gz",
        )

        # Test for malformed DE file. deseq2_bad_output.tab.gz file has a
        # mangled last value in last column where we replaced dot for a comma.
        # When importing to pandas DataFrame, this is no longer a float and
        # should throw an error.
        inputs["src"] = "./diff_exp/input/deseq2_bad_output.tab.gz"
        diff_bad = self.run_process("upload-diffexp", inputs, Data.STATUS_ERROR)
        error_msg = [
            "Column padj is not numeric. Please make sure that the input file has valid numeric values (i.e. "
            "periods for decimal places)."
        ]
        self.assertEqual(diff_bad.process_error, error_msg)

    @tag_process("upload-diffexp")
    def test_upload_de_check_field_type(self):
        inputs = {
            "src": "./diff_exp/input/diff_exp_check_geneid_type.tab.gz",
            "source": "DICTYBASE",
            "gene_id": "index",
            "logfc": "log2FoldChange",
            "fdr": "padj",
            "pvalue": "pvalue",
            "stat": "stat",
            "species": "Dictyostelium discoideum",
            "build": "dd-05-2009",
            "feature_type": "gene",
        }
        diff_exp = self.run_process("upload-diffexp", inputs)
        saved_json, test_json = self.get_json(
            "./diff_exp/output/diff_exp_check_types.json.gz", diff_exp.output["de_json"]
        )
        self.assertEqual(test_json, saved_json)
        all(self.assertIsInstance(gene, str) for gene in test_json["gene_id"])

    @tag_process("upload-bed")
    def test_upload_bed(self):
        inputs = {"src": "bad.bed", "species": "Homo sapiens", "build": "hg19"}
        bed = self.run_process("upload-bed", inputs, Data.STATUS_ERROR)

        inputs = {"src": "good.bed", "species": "Homo sapiens", "build": "hg19"}
        bed = self.run_process("upload-bed", inputs)
        self.assertFile(bed, "bed", "good.bed")
        self.assertFile(bed, "tbi_jbrowse", "good.bed.gz.tbi")

    @tag_process("upload-geneset")
    def test_upload_geneset(self):
        inputs = {"src": "geneset.tab.gz", "source": "UCSC", "species": "Homo sapiens"}
        geneset = self.run_process("upload-geneset", inputs)

        self.assertFile(geneset, "geneset", "geneset_out.tab.gz", compression="gzip")
        self.assertFields(geneset, "source", "UCSC")
        self.assertFields(geneset, "species", "Homo sapiens")
        self.assertJSON(geneset, geneset.output["geneset_json"], "", "geneset.json.gz")

    @tag_process("create-geneset")
    def test_create_geneset(self):
        inputs = {
            "genes": ["ABC", "DEF", "GHI"],
            "source": "UCSC",
            "species": "Homo sapiens",
        }
        geneset = self.run_process("create-geneset", inputs)

        self.assertFile(geneset, "geneset", "geneset_2.tab.gz", compression="gzip")
        self.assertJSON(
            geneset, geneset.output["geneset_json"], "", "geneset_2.json.gz"
        )
        self.assertFields(geneset, "source", "UCSC")
        self.assertFields(geneset, "species", "Homo sapiens")

        inputs = {
            "genes": ["1", "3", "3", "2"],
            "source": "NCBI",
            "species": "Homo sapiens",
        }
        geneset_2 = self.run_process("create-geneset", inputs)

        self.assertFile(geneset_2, "geneset", "geneset_3.tab.gz", compression="gzip")
        self.assertJSON(
            geneset_2, geneset_2.output["geneset_json"], "", "geneset_3.json.gz"
        )
        self.assertEqual(geneset_2.process_warning[0], "Removed duplicated genes.")

    @tag_process("create-geneset-venn")
    def test_create_venn(self):
        inputs = {
            "genes": ["ABC", "GHI", "DEF"],
            "source": "UCSC",
            "venn": "venn.json.gz",
            "species": "Homo sapiens",
        }
        venn = self.run_process("create-geneset-venn", inputs)

        self.assertFile(venn, "geneset", "geneset_venn.tab.gz", compression="gzip")
        self.assertJSON(venn, venn.output["geneset_json"], "", "geneset_venn.json.gz")
        self.assertJSON(venn, venn.output["venn"], "", "venn.json.gz")
        self.assertFields(venn, "source", "UCSC")
        self.assertFields(venn, "species", "Homo sapiens")

    @tag_process("upload-fastq-single")
    def test_upload_reformating_single(self):
        inputs = {"src": ["old_encoding.fastq.gz"]}
        reads = self.run_process("upload-fastq-single", inputs)
        self.assertFiles(
            reads, "fastq", ["old_encoding_transformed.fastq.gz"], compression="gzip"
        )

    @tag_process("upload-fastq-paired")
    def test_upload_reformating_paired(self):
        inputs = {
            "src1": ["old_encoding.fastq.gz", "old_encoding1.fastq.gz"],
            "src2": ["old_encoding_R2.fastq.gz", "old_encoding1_R2.fastq.gz"],
        }
        reads = self.run_process("upload-fastq-paired", inputs)
        self.assertFiles(
            reads,
            "fastq",
            ["old_encoding_transformed.fastq.gz", "old_encoding1_transformed.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            reads,
            "fastq2",
            [
                "old_encoding_transformed_R2.fastq.gz",
                "old_encoding1_transformed_R2.fastq.gz",
            ],
            compression="gzip",
        )

    @tag_process("upload-master-file")
    def test_upload_master_file(self):
        inputs = {"src": "56G_masterfile_corrupted.txt", "panel_name": "56G panel, v2"}
        master_file = self.run_process("upload-master-file", inputs, Data.STATUS_ERROR)

        # Check for non-unique amplicon names
        inputs["src"] = "56G masterfile_dup_amplicon_names.txt.gz"
        master_file = self.run_process("upload-master-file", inputs, Data.STATUS_ERROR)

        # Wrong file suffix
        inputs["src"] = "amplicon_master_file_merged.bed"
        master_file = self.run_process("upload-master-file", inputs, Data.STATUS_ERROR)

        # Check if primer sequences are allowed also in lowercase
        inputs["src"] = "56G masterfile_lowercase_bases.txt.gz"
        self.run_process("upload-master-file", inputs)

        inputs["src"] = "56G_masterfile_170113.txt.gz"
        master_file = self.run_process("upload-master-file", inputs)

        self.assertFile(master_file, "bedfile", "amplicon_master_file_merged.bed")
        self.assertFile(
            master_file, "nomergebed", "amplicon_master_file_nomergebed.bed"
        )
        self.assertFile(
            master_file, "olapfreebed", "amplicon_master_file_olapfreebed.bed"
        )
        self.assertFile(master_file, "primers", "amplicon_primers.bed")
        self.assertFields(master_file, "panel_name", "56G panel, v2")

    @tag_process("upload-etc")
    def test_upload_etc(self):
        inputs = {"src": "etc_upload_input.xls"}
        etc = self.run_process("upload-etc", inputs)

        self.assertFile(etc, "etcfile", "test_etc.json.gz")

    @tag_process("upload-fasta-nucl")
    def test_upload_nucl_seq(self):
        seq = self.run_process(
            "upload-fasta-nucl",
            {
                "src": os.path.join("nucl_seq", "input", "genome.fasta.gz"),
                "species": "Dictyostelium discoideum",
                "build": "dicty_2.7",
            },
        )
        self.assertFile(
            seq,
            "fastagz",
            os.path.join("nucl_seq", "output", "genome.fasta.gz"),
            compression="gzip",
        )
        self.assertFile(
            seq, "fasta", os.path.join("nucl_seq", "output", "genome.fasta")
        )
        self.assertFile(
            seq, "fai", os.path.join("nucl_seq", "output", "genome.fasta.fai")
        )
        self.assertFile(
            seq, "fasta_dict", os.path.join("nucl_seq", "output", "genome.dict")
        )
        self.assertFields(seq, "species", "Dictyostelium discoideum")
        self.assertFields(seq, "build", "dicty_2.7")
        self.assertFields(seq, "num_seqs", 1)

        empty_input = {
            "src": os.path.join("nucl_seq", "input", "empty.fasta"),
            "species": "Dictyostelium discoideum",
            "build": "dicty_2.7",
        }
        empty = self.run_process("upload-fasta-nucl", empty_input, Data.STATUS_ERROR)
        error_msg = ["The uploaded .FASTA file empty.fasta contains no sequence data."]
        self.assertEqual(empty.process_error, error_msg)

        malformed_input = {
            "src": os.path.join("nucl_seq", "input", "malformed.fasta"),
            "species": "Dictyostelium discoideum",
            "build": "dicty_2.7",
        }
        malformed = self.run_process(
            "upload-fasta-nucl", malformed_input, Data.STATUS_ERROR
        )
        error_msg = [
            "Format error in the uploaded file malformed.fasta. Error in FASTA file at "
            "line 1: Expected '>' at beginning of record, but got 'foo'."
        ]
        self.assertEqual(malformed.process_error, error_msg)

        incomplete_input = {
            "src": os.path.join("nucl_seq", "input", "incomplete.fasta"),
            "species": "Dictyostelium discoideum",
            "build": "dicty_2.7",
        }
        incomplete = self.run_process(
            "upload-fasta-nucl", incomplete_input, Data.STATUS_ERROR
        )
        error_msg = [
            "The uploaded .FASTA file incomplete.fasta contains no sequence data."
        ]
        self.assertEqual(incomplete.process_error, error_msg)

    @tag_process("upload-variants-vcf")
    def test_upload_vcf(self):
        vcf = self.run_process(
            "upload-variants-vcf",
            {"src": "igv_human.lf.vcf", "species": "Homo sapiens", "build": "b37",},
        )
        self.assertFile(vcf, "vcf", "igv_human.lf.vcf.gz", compression="gzip")
        self.assertFileExists(vcf, "tbi")
        self.assertFields(vcf, "species", "Homo sapiens")
        self.assertFields(vcf, "build", "b37")

    @tag_process("upload-gff3")
    def test_upload_gff3(self):
        inputs = {
            "src": "PGSC upload.gff3",
            "build": "ST",
            "source": "PGSC",
            "species": "Solanum tuberosum",
        }
        upload_gff3 = self.run_process("upload-gff3", inputs)
        del upload_gff3.output["annot_sorted_track_jbrowse"][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFile(upload_gff3, "annot_sorted", "PGSC upload_sorted.gff3")
        self.assertFields(
            upload_gff3,
            "annot_sorted_idx_igv",
            {"file": "PGSC upload_sorted.gff3.idx", "total_size": 126},
        )
        self.assertFields(
            upload_gff3,
            "annot_sorted_track_jbrowse",
            {"refs": ["tracks/annotation", "seq", "names"], "file": "trackList.json"},
        )
        self.assertFields(upload_gff3, "species", "Solanum tuberosum")
        self.assertFields(upload_gff3, "build", "ST")

    @tag_process("upload-gtf")
    def test_upload_gtf(self):
        inputs = {
            "src": "Hs GRCh38_86 upload.gtf",
            "build": "hg19",
            "source": "ENSEMBL",
            "species": "Homo Sapiens",
        }
        upload_gtf = self.run_process("upload-gtf", inputs)
        del upload_gtf.output["annot_sorted_track_jbrowse"][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFile(upload_gtf, "annot_sorted", "Hs GRCh38_86 upload_sorted.gtf")
        self.assertFields(
            upload_gtf,
            "annot_sorted_idx_igv",
            {"file": "Hs GRCh38_86 upload_sorted.gtf.idx", "total_size": 116},
        )
        self.assertFields(
            upload_gtf,
            "annot_sorted_track_jbrowse",
            {"refs": ["tracks/annotation", "seq", "names"], "file": "trackList.json"},
        )
        self.assertFields(upload_gtf, "species", "Homo Sapiens")
        self.assertFields(upload_gtf, "build", "hg19")

    @tag_process("upload-sc-10x")
    def test_upload_sc_reads(self):
        inputs = {
            "barcodes": ["10x_S1_L001_R1_001.fastq.gz", "10x_S1_L002_R1_001.fastq.gz"],
            "reads": ["10x_S1_L001_R2_001.fastq.gz"],
        }
        wrong_mates = self.run_process("upload-sc-10x", inputs, Data.STATUS_ERROR)
        error_msg = ["The number of reads and barcodes fastqs must be the same."]
        self.assertEqual(wrong_mates.process_error, error_msg)

        inputs = {
            "barcodes": ["10x_S1_L001_R1_001.fastq.gz", "10x_S1_L002_R1_001.fastq.gz"],
            "reads": ["10x_S1_L001_R2_001.fastq.gz", "10x_S1_L002_R2_001.fastq.gz"],
        }
        reads = self.run_process("upload-sc-10x", inputs)

        self.assertFiles(
            reads,
            "barcodes",
            ["10x_S1_L001_R1_001.fastq.gz", "10x_S1_L002_R1_001.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            reads,
            "reads",
            ["10x_S1_L001_R2_001.fastq.gz", "10x_S1_L002_R2_001.fastq.gz"],
            compression="gzip",
        )

    @tag_process("upload-bedpe")
    def test_upload_bedpe(self):
        species = "Homo sapiens"
        build = "fake_genome_RSEM"
        in_file = "./annotation_bedpe/input/alk.bedpe"

        inputs_bedpe = {"src": in_file, "species": species, "build": build}

        bedpe = self.run_process("upload-bedpe", inputs_bedpe)

        self.assertFile(bedpe, "bedpe", in_file)
        self.assertFields(bedpe, "species", species)
        self.assertFields(bedpe, "build", build)
