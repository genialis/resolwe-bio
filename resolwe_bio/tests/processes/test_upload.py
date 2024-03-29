import os
from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase, skipUnlessLargeFiles


class UploadProcessorTestCase(KBBioProcessTestCase):
    @tag_process("upload-bam", "upload-bam-indexed")
    def test_bam_upload(self):
        base_input = Path("test_bam")
        input_folder = base_input / "input"
        output_folder = base_input / "output"

        inputs = {
            "src": input_folder / "alignment_name_sorted.bam",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        upload_bam = self.run_process("upload-bam", inputs)
        self.assertFile(
            obj=upload_bam,
            field_path="bam",
            fn=output_folder / "alignment_position_sorted.bam",
        )
        self.assertFile(
            obj=upload_bam,
            field_path="bai",
            fn=output_folder / "alignment_bam_upload_index.bai",
        )
        self.assertFile(
            obj=upload_bam,
            field_path="stats",
            fn=output_folder / "alignment_bam_upload_stats.txt",
        )
        self.assertFields(obj=upload_bam, path="species", value="Homo sapiens")
        self.assertFields(obj=upload_bam, path="build", value="hg19")

        inputs = {
            "src": output_folder / "alignment_position_sorted.bam",
            "src2": input_folder / "alignment_bam_upload_index.bam.bai",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        upload_bam = self.run_process("upload-bam-indexed", inputs, Data.STATUS_ERROR)
        self.assertEqual(
            upload_bam.process_error[0],
            "BAI should have the same name as BAM with .bai extension",
        )

        inputs = {
            "src": output_folder / "alignment_position_sorted.bam",
            "src2": output_folder / "alignment_position_sorted.bam.bai",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        upload_bam = self.run_process("upload-bam-indexed", inputs)
        self.assertFile(
            obj=upload_bam,
            field_path="bam",
            fn=output_folder / "alignment_position_sorted.bam",
        )
        self.assertFile(
            obj=upload_bam,
            field_path="bai",
            fn=output_folder / "alignment_position_sorted.bam.bai",
        )
        self.assertFile(
            obj=upload_bam,
            field_path="stats",
            fn=output_folder / "alignment_bam_upload_stats.txt",
        )
        self.assertFields(obj=upload_bam, path="species", value="Homo sapiens")
        self.assertFields(obj=upload_bam, path="build", value="hg19")

    @with_resolwe_host
    @tag_process("upload-expression")
    def test_upload_expression(self):
        input_folder = Path("test_upload_expression") / "input"
        output_folder = Path("test_upload_expression") / "output"

        inputs = {
            "exp_type": "TPM",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        self.run_process("upload-expression", inputs, Data.STATUS_ERROR)

        inputs = {
            "exp": input_folder / "exp_1_tpm.tab.gz",
            "rc": input_folder / "exp_1_rc.tab.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        self.run_process("upload-expression", inputs, Data.STATUS_ERROR)

        inputs = {
            "rc": input_folder / "exp_1_rc.tab.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_3 = self.run_process("upload-expression", inputs)
        self.assertFile(exp_3, "rc", output_folder / "exp_1_rc.tab.gz")
        self.assertFile(exp_3, "exp", output_folder / "exp_1_rc.tab.gz")
        self.assertFile(
            exp_3,
            "exp_set",
            output_folder / "exp_1_rc_expressions3.txt.gz",
            compression="gzip",
        )
        self.assertJSON(
            exp_3, exp_3.output["exp_json"], "", output_folder / "exp_1.json.gz"
        )
        self.assertFields(exp_3, "species", "Homo sapiens")
        self.assertFields(exp_3, "build", "hg19")
        self.assertFields(exp_3, "feature_type", "gene")

        inputs = {
            "exp": input_folder / "exp_1_tpm.tab.gz",
            "exp_type": "TPM",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_4 = self.run_process("upload-expression", inputs)
        self.assertFile(exp_4, "exp", output_folder / "exp_1_tpm.tab.gz")
        self.assertFile(
            exp_4,
            "exp_set",
            output_folder / "exp_1_tpm_expressions.txt.gz",
            compression="gzip",
        )

        inputs = {
            "rc": input_folder / "exp_1_rc.tab.gz",
            "exp": input_folder / "exp_1_tpm.tab.gz",
            "exp_type": "TPM",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_5 = self.run_process("upload-expression", inputs)
        self.assertFields(exp_5, "exp_type", "TPM")
        self.assertFile(exp_5, "exp", output_folder / "exp_1_tpm.tab.gz")
        self.assertFile(exp_5, "rc", output_folder / "exp_1_rc.tab.gz")
        self.assertFile(
            exp_5,
            "exp_set",
            output_folder / "exp_1_rc_expressions5.txt.gz",
            compression="gzip",
        )
        self.assertJSON(
            exp_5, exp_5.output["exp_json"], "", output_folder / "exp_1_norm.json.gz"
        )
        self.assertJSON(
            exp_5,
            exp_5.output["exp_set_json"],
            "",
            output_folder / "upload_exp_norm_set.json.gz",
        )

        inputs = {
            "rc": input_folder / "exp_mac_line_ending.txt.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_6 = self.run_process("upload-expression", inputs)
        self.assertJSON(
            exp_6, exp_6.output["exp_json"], "", output_folder / "exp.json.gz"
        )
        self.assertFile(exp_6, "rc", output_folder / "exp_mac_line_ending.tab.gz")
        self.assertFile(
            exp_6,
            "exp_set",
            output_folder / "exp_mac_line_ending_expressions.txt.gz",
            compression="gzip",
        )

        inputs = {
            "rc": input_folder / "exp_unix_line_ending.txt.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_7 = self.run_process("upload-expression", inputs)
        self.assertJSON(
            exp_7, exp_7.output["exp_json"], "", output_folder / "exp.json.gz"
        )

        inputs = {
            "rc": input_folder / "exp_windows_line_ending.txt.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_8 = self.run_process("upload-expression", inputs)
        self.assertJSON(
            exp_8, exp_8.output["exp_json"], "", output_folder / "exp.json.gz"
        )

        # Check handling of numerical feature_ids in expression file
        inputs = {
            "rc": input_folder / "mm_ncbi_exp.tab.gz",
            "exp_name": "Expression",
            "source": "NCBI",
            "species": "Mus musculus",
            "build": "hg19",
        }
        exp_9 = self.run_process("upload-expression", inputs)
        self.assertFile(
            exp_9,
            "exp_set",
            output_folder / "mm_ncbi_exp_set.txt.gz",
            compression="gzip",
        )

        inputs = {
            "exp": input_folder / "exp_1_tpm.tsv.gz",
            "exp_type": "TPM",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_10 = self.run_process("upload-expression", inputs)
        self.assertFile(exp_10, "exp", output_folder / "exp_1_tpm.tab.gz")
        self.assertFile(
            exp_10,
            "exp_set",
            output_folder / "exp_1_tpm_expressions.txt.gz",
            compression="gzip",
        )

        # Test files with wrong extension
        inputs = {
            "rc": input_folder / "exp_1_rc.xlsx.gz",
            "exp": input_folder / "exp_1_tpm.tab.gz",
            "exp_type": "TPM",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        self.run_process("upload-expression", inputs, Data.STATUS_ERROR)

        inputs = {
            "exp": input_folder / "exp_1_tpm.gz",
            "exp_type": "TPM",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        self.run_process("upload-expression", inputs, Data.STATUS_ERROR)

        inputs = {
            "rc": input_folder / "exp.1_rc.tab.gz",
            "exp_name": "Expression",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        exp_13 = self.run_process("upload-expression", inputs)
        self.assertFields(exp_13, "rc", {"file": "exp.1_rc.tab.gz", "total_size": 99})
        self.assertFields(exp_13, "exp", {"file": "exp.1_rc.tab.gz", "total_size": 99})

    @with_resolwe_host
    @tag_process("upload-cxb", "upload-expression-cuffnorm")
    def test_upload_cuffquant_expr(self):
        input_folder = Path("test_upload_expression") / "input"
        output_folder = Path("test_upload_expression") / "output"
        inputs = {
            "src": "cuffquant 1.cxb",
            "source": "UCSC",
            "species": "Homo sapiens",
            "build": "hg19",
        }
        cxb = self.run_process("upload-cxb", inputs)

        inputs = {"exp": input_folder / "cuffquant_exp.tab", "cxb": cxb.id}
        exp = self.run_process("upload-expression-cuffnorm", inputs)
        self.assertFields(exp, "feature_type", "gene")
        self.assertFile(
            exp,
            "exp_set",
            output_folder / "cuffquant_exp_expressions.txt.gz",
            compression="gzip",
        )

    @tag_process("upload-fastq-paired")
    def test_upload_paired_end_reads(self):
        input_folder = Path("test_fastq_upload") / "input"
        output_folder = Path("test_fastq_upload") / "output"
        inputs = {
            "src1": [
                input_folder / "mate1_reordered.fastq.gz",
                input_folder / "mate1_diff_num_reads.fastq.gz",
            ],
            "src2": [input_folder / "mate2_reordered.fastq.gz"],
        }
        wrong_mates = self.run_process("upload-fastq-paired", inputs, Data.STATUS_ERROR)
        error_msg = [
            "The number of mate-pair files in split-lane samples must match. 2 and 1 "
            "input files were given for the -fq and -fq2 inputs, respectively."
        ]
        self.assertEqual(wrong_mates.process_error, error_msg)

        inputs = {
            "src1": [
                input_folder / "mate1_reordered.fastq.gz",
                input_folder / "mate1_reordered.fastq.gz",
            ],
            "src2": [input_folder / "mate2_reordered.fastq.gz"],
        }
        wrong_mates2 = self.run_process(
            "upload-fastq-paired", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            "Non-unique input file names detected: ['mate1_reordered.fastq.gz']."
        ]
        self.assertEqual(wrong_mates2.process_error, error_msg)

        inputs = {
            "src1": [input_folder / "mate1_diff_num_reads.fastq.gz"],
            "src2": [input_folder / "mate2_diff_num_reads.fastq.gz"],
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
            "src1": [input_folder / "mate1_reordered.fastq.gz"],
            "src2": [input_folder / "mate2_reordered.fastq.gz"],
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

        inputs = {
            "src1": [input_folder / "mate1.fastq.gz"],
            "src2": [input_folder / "mate2.fastq.gz"],
        }
        self.run_process("upload-fastq-paired", inputs, Data.STATUS_ERROR)

        inputs = {
            "src1": [input_folder / "rRNA_forw.fastq.gz"],
            "src2": [input_folder / "rRNA_rew.fastq.gz"],
        }
        reads = self.run_process("upload-fastq-paired", inputs)
        self.assertFiles(
            reads, "fastq", [output_folder / "rRNA_forw.fastq.gz"], compression="gzip"
        )
        self.assertFiles(
            reads, "fastq2", [output_folder / "rRNA_rew.fastq.gz"], compression="gzip"
        )
        del reads.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        self.assertFields(
            reads,
            "fastqc_url",
            [
                {
                    "file": "fastqc/rRNA_forw_fastqc/fastqc_report.html",
                    "refs": ["fastqc/rRNA_forw_fastqc"],
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
                "src1": [
                    input_folder / "old_encoding.fastq.gz",
                    input_folder / "old_encoding1.fastq.gz",
                ],
                "src2": [
                    input_folder / "old_encoding_R2.fastq.gz",
                    input_folder / "old_encoding1_R2.fastq.gz",
                ],
                "merge_lanes": True,
            },
        )
        self.assertFiles(
            merged_lanes,
            "fastq",
            [output_folder / "paired_end_merged_lanes_mate1.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            merged_lanes,
            "fastq2",
            [output_folder / "paired_end_merged_lanes_mate2.fastq.gz"],
            compression="gzip",
        )

        inputs = {
            "src1": [input_folder / "rRNA_forw.fastq"],
            "src2": [input_folder / "rRNA_rew.fastq"],
        }
        reads = self.run_process("upload-fastq-paired", inputs)
        self.assertFiles(
            reads, "fastq", [output_folder / "rRNA_forw.fastq.gz"], compression="gzip"
        )
        self.assertFiles(
            reads, "fastq2", [output_folder / "rRNA_rew.fastq.gz"], compression="gzip"
        )
        del reads.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        self.assertFields(
            reads,
            "fastqc_url",
            [
                {
                    "file": "fastqc/rRNA_forw_fastqc/fastqc_report.html",
                    "refs": ["fastqc/rRNA_forw_fastqc"],
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
        inputs = {
            "src1": [input_folder / "genome.fasta.gz"],
            "src2": [input_folder / "rRNA_rew.fastq"],
        }
        reads = self.run_process("upload-fastq-paired", inputs, Data.STATUS_ERROR)

    @tag_process("upload-fastq-single")
    def test_upload_single_end_reads(self):
        input_folder = Path("test_fastq_upload") / "input"
        output_folder = Path("test_fastq_upload") / "output"
        empty_input = self.run_process(
            "upload-fastq-single",
            {"src": [input_folder / "empty.fastq.gz"]},
            Data.STATUS_ERROR,
        )
        error_msg = ["Input file empty.fastq.gz contains no read sequences."]
        self.assertEqual(empty_input.process_error, error_msg)

        garbage_input = self.run_process(
            "upload-fastq-single",
            {"src": [input_folder / "garbage.fastq.gz"]},
            Data.STATUS_ERROR,
        )

        error_msg = [
            "Error in file garbage.fastq.gz. Error in FASTQ file at line 2: Premature "
            "end of file encountered. The incomplete final record was: 'Some random "
            "content\\n'"
        ]
        self.assertEqual(garbage_input.process_error, error_msg)

        garbage_input_2 = self.run_process(
            "upload-fastq-single",
            {"src": [input_folder / "garbage2.fastq.gz"]},
            Data.STATUS_ERROR,
        )
        error_msg = [
            "Error in file garbage2.fastq.gz. Error in FASTQ file at line 3: Sequence descriptions "
            "don't match ('Some random content' != '+++').\nThe second sequence description must "
            "be either empty or equal to the first description."
        ]
        self.assertEqual(garbage_input_2.process_error, error_msg)

        missing_qual = self.run_process(
            "upload-fastq-single",
            {"src": [input_folder / "missing_qual.fastq.gz"]},
            Data.STATUS_ERROR,
        )
        error_msg = [
            "Error in file missing_qual.fastq.gz. Error in FASTQ file at line 16: "
            "Premature end of file encountered. The incomplete final record was: "
            "'@read4/1\\nGACAGGCCGTTTGAATGTTGACGGGATGTT\\n+\\n'"
        ]
        self.assertEqual(missing_qual.process_error, error_msg)

        inputs = {"src": [input_folder / "mate1.fastq.gz"]}
        self.run_process("upload-fastq-single", inputs, Data.STATUS_ERROR)

        inputs = {
            "src": [
                input_folder / "rRNA_forw.fastq.gz",
                input_folder / "rRNA_rew.fastq.gz",
            ]
        }
        reads = self.run_process("upload-fastq-single", inputs)

        self.assertFiles(
            reads,
            "fastq",
            [
                output_folder / "rRNA_forw_single.fastq.gz",
                output_folder / "rRNA_rew.fastq.gz",
            ],
            compression="gzip",
        )
        del reads.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        del reads.output["fastqc_url"][1]["total_size"]  # Non-deterministic output.
        self.assertFields(
            reads,
            "fastqc_url",
            [
                {
                    "file": "fastqc/rRNA_forw_fastqc/fastqc_report.html",
                    "refs": ["fastqc/rRNA_forw_fastqc"],
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
            {
                "src": [
                    input_folder / "rRNA_forw.fastq.gz",
                    input_folder / "rRNA_rew.fastq.gz",
                ],
                "merge_lanes": True,
            },
        )
        self.assertFiles(
            merged_lanes,
            "fastq",
            [output_folder / "merged_single_end_reads.fastq.gz"],
            compression="gzip",
        )

        inputs = {
            "src": [
                input_folder / "rRNA_forw.fq.gz",
                input_folder / "rRNA_rew.fq.gz",
            ]
        }
        reads = self.run_process("upload-fastq-single", inputs)

        self.assertFiles(
            reads,
            "fastq",
            [
                output_folder / "rRNA_forw.fastq.gz",
                output_folder / "rRNA_rew.fastq.gz",
            ],
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
        self.assertEqual(geneset.descriptor_schema.slug, "geneset")

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
        self.assertEqual(geneset.descriptor_schema.slug, "geneset")

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
        self.assertEqual(geneset_2.descriptor_schema.slug, "geneset")

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
        self.assertEqual(venn.descriptor_schema.slug, "geneset")

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

    @tag_process("upload-etc")
    def test_upload_etc(self):
        inputs = {"src": "etc_upload_input.xls"}
        etc = self.run_process("upload-etc", inputs)

        self.assertFile(etc, "etcfile", "test_etc.json.gz")

    @tag_process("upload-fasta-nucl")
    def test_upload_nucl_seq(self):
        self.run_process(
            "upload-fasta-nucl",
            {
                "src": os.path.join("nucl_seq", "input", "chrX_1_30000.fasta.gz"),
                "species": "Homo sapiens",
                "build": "test",
            },
        )

        self.run_process(
            "upload-fasta-nucl",
            {
                "src": os.path.join("nucl_seq", "input", "chrX_1_30000.fasta"),
                "species": "Homo sapiens",
                "build": "test",
            },
        )

        self.run_process(
            "upload-fasta-nucl",
            {
                "src": os.path.join("nucl_seq", "input", "chrX_1_30000.fa.gz"),
                "species": "Homo sapiens",
                "build": "test",
            },
        )

        wrong_extension = self.run_process(
            "upload-fasta-nucl",
            {
                "src": os.path.join("nucl_seq", "input", "chrX_1_30000.fastq.gz"),
                "species": "Homo sapiens",
                "build": "test",
            },
            Data.STATUS_ERROR,
        )
        error_msg = [
            "The imported file has unsupported file name extension. "
            "The supported extensions are ('.fa', '.fasta')."
        ]
        self.assertEqual(wrong_extension.process_error, error_msg)

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
            {
                "src": "igv_human.lf.vcf",
                "species": "Homo sapiens",
                "build": "b37",
            },
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
        del upload_gff3.output["annot_sorted_idx_igv"][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFields(
            upload_gff3,
            "annot_sorted_idx_igv",
            {"file": "PGSC upload_sorted.gff3.idx"},
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
        del upload_gtf.output["annot_sorted_idx_igv"][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFields(
            upload_gtf,
            "annot_sorted_idx_igv",
            {"file": "Hs GRCh38_86 upload_sorted.gtf.idx"},
        )
        self.assertFields(
            upload_gtf,
            "annot_sorted_track_jbrowse",
            {"refs": ["tracks/annotation", "seq", "names"], "file": "trackList.json"},
        )
        self.assertFields(upload_gtf, "species", "Homo Sapiens")
        self.assertFields(upload_gtf, "build", "hg19")

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

    @tag_process("upload-proteomics-sample")
    def test_upload_proteomics_sample(self):
        base = Path("proteomics")
        inputs = base / "input"
        outputs = base / "output"

        prot_data = self.run_process(
            "upload-proteomics-sample",
            {
                "src": str(inputs / "single_sample.txt"),
                "species": "Homo sapiens",
            },
        )
        self.assertFile(prot_data, "table", str(outputs / "single_sample.txt"))
        self.assertFields(prot_data, "species", "Homo sapiens")
        self.assertFields(prot_data, "source", "UniProtKB")

    @tag_process("upload-proteomics-sample-set")
    def test_upload_proteomics_sample_set(self):
        base = Path("proteomics")
        inputs = base / "input"
        outputs = base / "output"

        prot_data = self.run_process(
            "upload-proteomics-sample-set",
            {
                "src": str(inputs / "sample_set.txt"),
                "species": "Homo sapiens",
            },
        )
        self.assertFile(prot_data, "table", str(outputs / "sample_set.txt"))
        self.assertFields(prot_data, "species", "Homo sapiens")
        self.assertFields(prot_data, "source", "UniProtKB")

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertEqual(
            Data.objects.filter(process__slug="upload-proteomics-sample").count(), 2
        )

    @skipUnlessLargeFiles(
        "6042316072_R03C01_Red.idat.gz", "6042316072_R03C01_Grn.idat.gz"
    )
    @tag_process("upload-idat")
    def test_upload_idat(self):
        large = Path("large")
        base_path = Path("methylation") / "inputs"
        inputs = {
            "red_channel": large / "6042316072_R03C01_Red.idat.gz",
            "green_channel": large / "6042316072_R03C01_Grn.idat.gz",
            "species": "Homo sapiens",
            "platform": "HM450",
        }

        idat = self.run_process("upload-idat", inputs)

        self.assertFile(
            idat,
            "red_channel",
            large / "6042316072_R03C01_Red.idat.gz",
            compression="gzip",
        )
        self.assertFile(
            idat,
            "green_channel",
            large / "6042316072_R03C01_Grn.idat.gz",
            compression="gzip",
        )
        self.assertFields(idat, "species", "Homo sapiens")
        self.assertFields(idat, "platform", "HM450")

        inputs.update({"species": "Mus musculus"})
        wrong_species_input = self.run_process("upload-idat", inputs, Data.STATUS_ERROR)
        error_msg = [
            "Platform type HM450 does not match the selected species Mus musculus."
        ]
        self.assertEqual(wrong_species_input.process_error, error_msg)

        inputs.update(
            {
                "red_channel": base_path / "wrong_sample_name_Blue.idat.gz",
                "species": "Homo sapiens",
            }
        )
        self.run_process("upload-idat", inputs, Data.STATUS_ERROR)

    @tag_process("upload-vep-cache")
    def test_upload_vep_cache(self):
        input_folder = Path("ensembl-vep") / "input"
        output_folder = Path("ensembl-vep") / "output"
        vep_cache = self.run_process(
            "upload-vep-cache",
            {
                "cache_file": input_folder / "cache_homo_sapiens_X.tar.gz",
                "species": "Homo sapiens",
                "build": "GRCh38",
                "release": "104",
            },
        )
        self.assertDir(
            vep_cache, "cache", output_folder / "cache_homo_sapiens_X.tar.gz"
        )
        self.assertFields(vep_cache, "species", "Homo sapiens")
        self.assertFields(vep_cache, "build", "GRCh38")
        self.assertFields(vep_cache, "release", "104")
