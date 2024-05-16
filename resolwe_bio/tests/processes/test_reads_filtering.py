from pathlib import Path

from django.test import override_settings

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_docker_executor

from resolwe_bio.utils.test import BioProcessTestCase


class ReadsFilteringProcessorTestCase(BioProcessTestCase):
    @tag_process("trimmomatic-single")
    def test_trimmomatic_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads()
            adapters = self.run_process(
                "upload-fasta-nucl",
                {
                    "src": Path("bbduk", "input", "bbduk_adapters.fasta"),
                    "species": "Other",
                    "build": "Illumina adapters",
                },
            )

        inputs = {
            "reads": reads.pk,
            "illuminaclip": {
                "adapters": adapters.pk,
                "seed_mismatches": 2,
                "simple_clip_threshold": 10,
            },
            "maxinfo": {
                "target_length": 10,
                "strictness": 0.6,
            },
            "slidingwindow": {
                "window_size": 4,
                "required_quality": 15,
            },
            "trim_bases": {
                "leading": 20,
                "trailing": 20,
                "crop": 40,
                "headcrop": 3,
            },
            "reads_filtering": {
                "minlen": 22,
                "average_quality": 10,
            },
        }
        filtered_reads = self.run_processor("trimmomatic-single", inputs)

        self.assertFiles(
            filtered_reads,
            "fastq",
            ["filtered_reads_trimmomatic_single.fastq.gz"],
            compression="gzip",
        )
        del filtered_reads.output["fastqc_url"][0][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFields(
            filtered_reads,
            "fastqc_url",
            [
                {
                    "file": "fastqc/reads_fastqc/fastqc_report.html",
                    "refs": ["fastqc/reads_fastqc"],
                }
            ],
        )

    @tag_process("trimmomatic-paired")
    def test_trimmomatic_paired(self):
        with self.preparation_stage():
            inputs = {"src1": ["rRNA_forw.fastq.gz"], "src2": ["rRNA_rew.fastq.gz"]}
            reads = self.run_processor("upload-fastq-paired", inputs)

        inputs = {"reads": reads.pk, "trim_bases": {"trailing": 3}}

        filtered_reads = self.run_processor("trimmomatic-paired", inputs)
        self.assertFiles(
            filtered_reads,
            "fastq",
            ["filtered_reads_trimmomatic_paired_fw.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            filtered_reads,
            "fastq2",
            ["filtered_reads_trimmomatic_paired_rw.fastq.gz"],
            compression="gzip",
        )
        del filtered_reads.output["fastqc_url"][0][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFields(
            filtered_reads,
            "fastqc_url",
            [
                {
                    "file": "fastqc/rRNA_forw_fastqc/fastqc_report.html",
                    "refs": ["fastqc/rRNA_forw_fastqc"],
                }
            ],
        )
        del filtered_reads.output["fastqc_url2"][0][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFields(
            filtered_reads,
            "fastqc_url2",
            [
                {
                    "file": "fastqc/rRNA_rew_fastqc/fastqc_report.html",
                    "refs": ["fastqc/rRNA_rew_fastqc"],
                }
            ],
        )

    @tag_process("cutadapt-single")
    def test_cutadapt_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads(
                ["cutadapt single.fastq.gz", "cutadapt_single1.fastq.gz"]
            )

            primers_up = self.prepare_ref_seq("5_prime_adapter.fasta.gz")
            primers_down = self.prepare_ref_seq("3_prime_adapter.fasta.gz")

        inputs = {
            "reads": reads.id,
            "adapters": {
                "polya_tail": 5,
                "down_primers_seq": ["AGCACCT"],
                "up_primers_seq": ["AGCTAAA"],
            },
            "modify_reads": {
                "nextseq_trim": 5,
            },
            "filtering": {
                "minlen": 10,
            },
        }

        cutadapt_single = self.run_process("cutadapt-single", inputs)

        self.assertFiles(
            cutadapt_single,
            "fastq",
            ["cutadapt_single_trimmed.fastq.gz"],
            compression="gzip",
        )

        inputs = {
            "reads": reads.id,
            "adapters": {
                "polya_tail": 5,
                "down_primers_file": primers_down.id,
                "up_primers_file": primers_up.id,
            },
            "filtering": {
                "minlen": 10,
            },
        }

        cutadapt_single = self.run_process("cutadapt-single", inputs)

        self.assertFiles(
            cutadapt_single,
            "fastq",
            ["cutadapt_single_trimmed.fastq.gz"],
            compression="gzip",
        )

    @tag_process("cutadapt-paired")
    def test_cutadapt_paired(self):
        with self.preparation_stage():
            reads = self.prepare_paired_reads(
                mate1=["cutadapt mate1.fastq.gz"], mate2=["cutadapt mate2.fastq.gz"]
            )

            primers_up = self.prepare_ref_seq("5_prime_adapter.fasta.gz")
            primers_down = self.prepare_ref_seq("3_prime_adapter.fasta.gz")

        inputs = {
            "reads": reads.id,
            "adapters": {
                "mate1_3prime_seq": ["AGCACCT"],
                "mate2_3prime_seq": ["AGCACCT"],
                "mate1_5prime_seq": ["AGCTAAA"],
                "mate2_5prime_seq": ["AGCTAAA"],
            },
            "filtering": {
                "minlen": 10,
            },
        }

        cutadapt_paired = self.run_process("cutadapt-paired", inputs)

        self.assertFiles(
            cutadapt_paired,
            "fastq",
            ["cutadapt_paired_forward_trimmed.fastq.gz"],
            compression="gzip",
        )

        self.assertFiles(
            cutadapt_paired,
            "fastq2",
            ["cutadapt_paired_reverse_trimmed.fastq.gz"],
            compression="gzip",
        )

        inputs = {
            "reads": reads.id,
            "adapters": {
                "mate1_3prime_file": primers_down.id,
                "mate2_3prime_file": primers_down.id,
                "mate1_5prime_file": primers_up.id,
                "mate2_5prime_file": primers_up.id,
            },
            "filtering": {
                "minlen": 10,
            },
        }

        cutadapt_paired = self.run_process("cutadapt-paired", inputs)

        self.assertFiles(
            cutadapt_paired,
            "fastq",
            ["cutadapt_paired_forward_trimmed.fastq.gz"],
            compression="gzip",
        )

        self.assertFiles(
            cutadapt_paired,
            "fastq2",
            ["cutadapt_paired_reverse_trimmed.fastq.gz"],
            compression="gzip",
        )

    @tag_process("cutadapt-3prime-single")
    def test_cutadapt_3prime_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads(
                ["cutadapt single.fastq.gz", "cutadapt_single1.fastq.gz"]
            )

        inputs = {
            "reads": reads.id,
            "options": {
                "nextseq_trim": 5,
                "min_len": 20,
                "min_overlap": 20,
                "times": 2,
            },
        }
        cutadapt_single = self.run_process("cutadapt-3prime-single", inputs)

        self.assertFiles(
            cutadapt_single,
            "fastq",
            ["cutadapt_3prime_single_trimmed.fastq.gz"],
            compression="gzip",
        )

    @tag_process("cutadapt-corall-single")
    def test_cutadapt_corall_single(self):
        with self.preparation_stage():
            reads = self.prepare_reads(["./corall/input/corall_single.fastq.gz"])

        cutadapt_single = self.run_process(
            "cutadapt-corall-single", {"reads": reads.id}
        )

        self.assertFiles(
            cutadapt_single,
            "fastq",
            ["./corall/output/single_trimmed.fastq.gz"],
            compression="gzip",
        )

    @tag_process("cutadapt-corall-paired")
    def test_cutadapt_corall_paired(self):
        with self.preparation_stage():
            reads_paired = self.prepare_paired_reads(
                mate1=["./corall/input/corall_mate1.fastq.gz"],
                mate2=["./corall/input/corall_mate2.fastq.gz"],
            )

        cutadapt_paired = self.run_process(
            "cutadapt-corall-paired", {"reads": reads_paired.id}
        )

        self.assertFiles(
            cutadapt_paired,
            "fastq",
            ["./corall/output/mate1_trimmed.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            cutadapt_paired,
            "fastq2",
            ["./corall/output/mate2_trimmed.fastq.gz"],
            compression="gzip",
        )

    @with_docker_executor
    @override_settings(FLOW_PROCESS_MAX_CORES=4)
    @tag_process("bbduk-single", "bbduk-paired")
    def test_bbduk(self):
        input_folder = Path("bbduk") / "input"
        output_folder = Path("bbduk") / "output"
        with self.preparation_stage():
            reads = self.prepare_reads(
                [
                    input_folder / "bbduk_test_reads.fastq.gz",
                    input_folder / "rRNA_forw.fastq.gz",
                    "fw reads.fastq.gz",
                ]
            )
            reads_single = self.prepare_reads(
                [input_folder / "bbduk_test_reads.fastq.gz"]
            )
            reads_paired = self.prepare_paired_reads(
                [input_folder / "rRNA_forw.fastq.gz"],
                [input_folder / "rRNA_rew.fastq.gz"],
            )
            ref_seq = self.prepare_ref_seq(
                fn=input_folder / "bbduk_adapters.fasta",
                species="Other",
                build="Custom",
            )
            barcodes = self.prepare_ref_seq(
                fn=input_folder / "barcodes.fasta",
                species="Other",
                build="Custom",
            )
            reads_paired_lanes = self.prepare_paired_reads(
                mate1=["fw reads.fastq.gz", "fw reads_2.fastq.gz"],
                mate2=["rw reads.fastq.gz", "rw reads_2.fastq.gz"],
            )

        inputs = {
            "reads": reads.id,
            "reference": {
                "sequences": [ref_seq.id],
            },
            "header_parsing": {
                "barcode_files": [barcodes.id],
            },
            "operations": {
                "quality_trim": "l",
            },
        }
        filtered_reads = self.run_process("bbduk-single", inputs)

        self.assertFiles(
            filtered_reads,
            "fastq",
            [
                output_folder / "bbduk_test_reads_preprocessed.fastq.gz",
                output_folder / "rRNA_forw_preprocessed.fastq.gz",
                output_folder / "fw_reads_preprocessed.fastq.gz",
            ],
            compression="gzip",
        )
        del filtered_reads.output["fastqc_url"][0][
            "total_size"
        ]  # Non-deterministic output.
        del filtered_reads.output["fastqc_url"][1]["total_size"]
        del filtered_reads.output["fastqc_url"][2]["total_size"]
        report = [
            {
                "file": "fastqc/bbduk_test_reads_preprocessed_fastqc/fastqc_report.html",
                "refs": [
                    "fastqc/bbduk_test_reads_preprocessed_fastqc",
                ],
            },
            {
                "file": "fastqc/rRNA_forw_preprocessed_fastqc/fastqc_report.html",
                "refs": ["fastqc/rRNA_forw_preprocessed_fastqc"],
            },
            {
                "file": "fastqc/fw reads_preprocessed_fastqc/fastqc_report.html",
                "refs": ["fastqc/fw reads_preprocessed_fastqc"],
            },
        ]
        self.assertFields(filtered_reads, "fastqc_url", report)

        inputs["reads"] = reads_single.id
        filtered_reads = self.run_process("bbduk-single", inputs)

        inputs = {
            "reads": reads_paired.id,
        }

        filtered_reads = self.run_process("bbduk-paired", inputs)

        self.assertFiles(
            filtered_reads,
            "fastq",
            [output_folder / "bbduk_fw_reads.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            filtered_reads,
            "fastq2",
            [output_folder / "bbduk_rv_reads.fastq.gz"],
            compression="gzip",
        )
        del filtered_reads.output["fastqc_url"][0][
            "total_size"
        ]  # Non-deterministic output.
        report = {
            "file": "fastqc/rRNA_forw_preprocessed_fastqc/fastqc_report.html",
            "refs": [
                "fastqc/rRNA_forw_preprocessed_fastqc",
            ],
        }
        self.assertFields(filtered_reads, "fastqc_url", [report])
        del filtered_reads.output["fastqc_url2"][0][
            "total_size"
        ]  # Non-deterministic output.
        report2 = {
            "file": "fastqc/rRNA_rew_preprocessed_fastqc/fastqc_report.html",
            "refs": [
                "fastqc/rRNA_rew_preprocessed_fastqc",
            ],
        }
        self.assertFields(filtered_reads, "fastqc_url2", [report2])

        inputs = {
            "reads": reads_paired_lanes.id,
        }

        filtered_reads = self.run_process("bbduk-paired", inputs)

        del filtered_reads.output["fastqc_url"][0][
            "total_size"
        ]  # Non-deterministic output.
        del filtered_reads.output["fastqc_url"][1][
            "total_size"
        ]  # Non-deterministic output.
        report = [
            {
                "file": "fastqc/fw reads_preprocessed_fastqc/fastqc_report.html",
                "refs": [
                    "fastqc/fw reads_preprocessed_fastqc",
                ],
            },
            {
                "file": "fastqc/fw reads_2_preprocessed_fastqc/fastqc_report.html",
                "refs": [
                    "fastqc/fw reads_2_preprocessed_fastqc",
                ],
            },
        ]
        self.assertFields(filtered_reads, "fastqc_url", report)

        del filtered_reads.output["fastqc_url2"][0][
            "total_size"
        ]  # Non-deterministic output.
        del filtered_reads.output["fastqc_url2"][1][
            "total_size"
        ]  # Non-deterministic output.
        report2 = [
            {
                "file": "fastqc/rw reads_preprocessed_fastqc/fastqc_report.html",
                "refs": [
                    "fastqc/rw reads_preprocessed_fastqc",
                ],
            },
            {
                "file": "fastqc/rw reads_2_preprocessed_fastqc/fastqc_report.html",
                "refs": [
                    "fastqc/rw reads_2_preprocessed_fastqc",
                ],
            },
        ]
        self.assertFields(filtered_reads, "fastqc_url2", report2)

    @tag_process("bamclipper")
    def test_bamclipper(self):
        species = "Homo sapiens"
        build = "fake_genome_RSEM"
        align_input = "./bamclipper/input/TP53.bam"

        with self.preparation_stage():
            bam = self.prepare_bam(fn=align_input, species=species, build=build)

            inputs_bedpe = {
                "src": "./bamclipper/input/TP53.bedpe",
                "species": species,
                "build": build,
            }
            bedpe = self.run_process("upload-bedpe", inputs_bedpe)

        # Test if bamclipper has been skipped.
        bc_skip_inputs = {"alignment": bam.id, "skip": True}
        skipped_bc = self.run_process("bamclipper", bc_skip_inputs)
        self.assertFile(skipped_bc, "bam", align_input)
        bc_data = Data.objects.last()
        self.assertEqual(bc_data.process_info, ["Skipping bamclipper step."])

        # Test bamclipper.
        inputs_bamclipper = {"alignment": bam.id, "bedpe": bedpe.id}
        clipped = self.run_process("bamclipper", inputs_bamclipper)

        self.assertFile(
            clipped, "stats", "./bamclipper/output/TP53.primerclipped.bam_stats.txt"
        )
        self.assertFields(clipped, "species", species)
        self.assertFields(clipped, "build", build)

    @tag_process("markduplicates")
    def test_markduplicates(self):
        species = "Homo sapiens"
        build = "custombuild"
        primerclipped = "./bamclipper/output/TP53.primerclipped.bam"

        with self.preparation_stage():
            bam = self.prepare_bam(fn=primerclipped, species=species, build=build)

        # Test if skipped. Input bam should always equal output bam.
        md_inputs = {"bam": bam.id, "skip": True}
        skipped_md = self.run_process("markduplicates", md_inputs)
        self.assertFile(skipped_md, "bam", primerclipped)

        # Test that removal of duplicates works.
        md_inputs = {
            "bam": bam.id,
            "remove_duplicates": True,
            "advanced": {
                "java_gc_threads": 3,
                "max_heap_size": 10,
            },
        }
        removed_md = self.run_process("markduplicates", md_inputs)

        def filter_startedon(line):
            return line.startswith(b"# Started on:") or line.startswith(
                b"# MarkDuplicates"
            )

        self.assertFileExists(removed_md, "bam")
        self.assertFileExists(removed_md, "bai")
        self.assertFile(
            removed_md,
            "stats",
            "./markduplicate/output/TP53.primerclipped.markduplicates.bam_stats.txt",
        )
        self.assertFile(
            removed_md,
            "metrics_file",
            "./markduplicate/output/TP53.primerclipped_metrics.txt",
            file_filter=filter_startedon,
        )
        self.assertFields(removed_md, "species", species)
        self.assertFields(removed_md, "build", build)

    @tag_process("bqsr")
    def test_bqsr(self):
        species = "Homo sapiens"
        build = "custom_build"
        with self.preparation_stage():
            input_genome = {
                # Based on b37 genome, chromosome 19 has been cut from beginning up to position 1207173.
                # This includes an exon of STK11. Cutting from the start of the chromosome was done so that
                # there is no need to shift any subsequent bed and vcf files.
                "src": "./bqsr/input/hs_b37_chr17_upto_TP53.fasta.gz",
                "species": species,
                "build": build,
            }
            input_bam = {
                "src": "./markduplicate/output/TP53.primerclipped.markduplicates.bam",
                "species": species,
                "build": build,
            }

            ks_dbsnp = []
            for i in ["./bqsr/input/dbsnp_TP53.vcf.gz"]:  # add more files if needed
                ks_dbsnp.append(
                    self.run_process(
                        "upload-variants-vcf",
                        {"src": i, "species": species, "build": build},
                    )
                )

            intervals = self.run_process(
                "upload-bed",
                {"src": "./bqsr/input/TP53.bed", "species": species, "build": build},
            )

            bam = self.run_process("upload-bam", input_bam)
            reference = self.run_process("upload-fasta-nucl", input_genome)

            bqsr_inputs = {
                "bam": bam.id,
                "reference": reference.id,
                "known_sites": [i.id for i in ks_dbsnp],
                "intervals": intervals.id,
            }
            bqsr = self.run_process("bqsr", bqsr_inputs)

            self.assertFileExists(bqsr, "bam")
            self.assertFileExists(bqsr, "bai")
            self.assertFile(
                bqsr,
                "stats",
                "./bqsr/output/TP53.primerclipped.markduplicates.bam_stats.txt",
            )
            self.assertFile(
                bqsr,
                "recal_table",
                "./bqsr/output/TP53.primerclipped.markduplicates_recalibration.table",
            )
            self.assertFields(bqsr, "species", species)
            self.assertFields(bqsr, "build", build)

            # Check if read groups has successfully been added.
            bqsr_inputs["read_group"] = "-LB=DAB;-PL=Illumina;-PU=barcode;-SM=sample1"
            bqsr_rg = self.run_process("bqsr", bqsr_inputs)

            self.assertFileExists(bqsr_rg, "bam")
            self.assertFileExists(bqsr_rg, "bai")

            bqsr_inputs["read_group"] = (
                "-LB=DAB;-PL=Illumina;-PU=barcode;-SM=sample1;-SM=sample2"
            )
            bqsr_dbltag = self.run_process("bqsr", bqsr_inputs, Data.STATUS_ERROR)
            self.assertEqual(
                bqsr_dbltag.process_error[0],
                "You have duplicate tags in read_group argument.",
            )

            bqsr_inputs["read_group"] = (
                "LB=DAB;-PL=Illumina;-PU=barcode;-SM=sample1"  # missing dash in LB
            )
            bqsr_tagerror = self.run_process("bqsr", bqsr_inputs, Data.STATUS_ERROR)
            self.assertEqual(
                bqsr_tagerror.process_error[0],
                "One or more read_group argument(s) improperly formatted.",
            )

            bqsr_inputs = {
                "bam": bam.id,
                "reference": reference.id,
                "known_sites": [i.id for i in ks_dbsnp],
                "intervals": intervals.id,
                "advanced": {
                    "use_original_qualities": True,
                    "java_gc_threads": 3,
                    "max_heap_size": 10,
                },
            }
            bqsr = self.run_process("bqsr", bqsr_inputs)
            self.assertFileExists(bqsr, "bam")
            self.assertFileExists(bqsr, "bai")
            self.assertFile(
                bqsr,
                "stats",
                "./bqsr/output/TP53.primerclipped.markduplicates.bam_stats.txt",
            )

    @tag_process("alignmentsieve")
    def test_alignmentsieve(self):
        species = "Homo sapiens"
        build = "custom_build"

        with self.preparation_stage():
            input_bam = {
                "src": "./markduplicate/output/TP53.primerclipped.markduplicates.bam",
                "species": species,
                "build": build,
            }
            bam = self.run_process("upload-bam", input_bam)

        params = {
            "alignment": bam.id,
            "max_fragment_length": 149,
        }

        # Produced .bam files are not deterministic and are not ideal
        # to run tests on. Ergo, we are testing for the stats file.
        max_filteredbam = self.run_process("alignmentsieve", params)
        self.assertFile(
            obj=max_filteredbam,
            field_path="stats",
            fn="./test_alignmentsieve/output/filtered_max_149_stats.txt",
        )

        params = {
            "alignment": bam.id,
            "min_fragment_length": 150,
        }
        min_filteredbam = self.run_process("alignmentsieve", params)
        self.assertFile(
            min_filteredbam,
            field_path="stats",
            fn="./test_alignmentsieve/output/filtered_min_150_stats.txt",
        )

    @tag_process("trimgalore-paired")
    def test_trimgalore_paired(self):
        with self.preparation_stage():
            reads = self.prepare_paired_reads(
                mate1=["cutadapt mate1.fastq.gz"], mate2=["cutadapt mate2.fastq.gz"]
            )
            adapters = self.prepare_ref_seq(
                "trimgalore/inputs/illumina_universal_adapters.fasta.gz"
            )

        inputs = {
            "reads": reads.id,
            "adapter_trim": {
                "universal_adapter": "--illumina",
            },
        }
        trimgalore = self.run_process("trimgalore-paired", inputs)

        self.assertFiles(
            trimgalore,
            "fastq",
            ["trimgalore/outputs/mate1_trimmed.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            trimgalore,
            "fastq2",
            ["trimgalore/outputs/mate2_trimmed.fastq.gz"],
            compression="gzip",
        )
        self.assertFileExists(trimgalore, "report")

        inputs = {
            "reads": reads.id,
            "adapter_trim": {
                "adapter": ["AGATCGGAAGAGC"],
                "adapter_2": ["AGATCGGAAGAGC"],
            },
        }
        trimgalore_adapters = self.run_process("trimgalore-paired", inputs)
        self.assertFiles(
            trimgalore_adapters,
            "fastq",
            ["trimgalore/outputs/mate1_trimmed.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            trimgalore_adapters,
            "fastq2",
            ["trimgalore/outputs/mate2_trimmed.fastq.gz"],
            compression="gzip",
        )

        inputs = {
            "reads": reads.id,
            "adapter_trim": {
                "adapter_file_1": adapters.id,
                "adapter_file_2": adapters.id,
            },
        }
        trimgalore_adapters = self.run_process("trimgalore-paired", inputs)
        self.assertFiles(
            trimgalore_adapters,
            "fastq",
            ["trimgalore/outputs/mate1_trimmed.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            trimgalore_adapters,
            "fastq2",
            ["trimgalore/outputs/mate2_trimmed.fastq.gz"],
            compression="gzip",
        )

    @tag_process("gatk-split-ncigar")
    def test_split_Ncigar_reads(self):
        input_folder = Path("splitNcigar_reads") / "input"
        output_folder = Path("splitNcigar_reads") / "output"
        with self.preparation_stage():
            species = "Homo sapiens"
            build = "GRCh38"

            bam = self.prepare_bam(
                fn=(input_folder / "chr1_500.bam"),
                species=species,
                build=build,
            )
            ref_seq = self.prepare_ref_seq(
                fn=(input_folder / "chr1_1-15000.fasta.gz"),
                species=species,
                build=build,
            )

        splitNcigar = self.run_process(
            "gatk-split-ncigar",
            {
                "bam": bam.id,
                "ref_seq": ref_seq.id,
                "advanced": {
                    "java_gc_threads": 3,
                    "max_heap_size": 10,
                },
            },
        )

        self.assertFileExists(splitNcigar, "bam")
        self.assertFileExists(splitNcigar, "bai")
        self.assertFile(
            splitNcigar,
            "stats",
            output_folder / "split_cigar_stats.txt",
        )
        self.assertFields(splitNcigar, "species", species)
        self.assertFields(splitNcigar, "build", build)

    @tag_process("xengsort-index", "xengsort-classify")
    def test_xengsort(self):
        def filter_variable_lines(line):
            """Filter variable lines."""
            if line.startswith(b"#"):
                return True
            elif b"time" in line:
                return True
            elif b"Time" in line:
                return True
            elif b"Computing weak k-mers" in line:
                return True
            elif b"hash files" in line:
                return True
            elif b"Done." in line:
                return True

        input_folder = Path("xengsort") / "input"
        output_folder = Path("xengsort") / "output"
        with self.preparation_stage():
            graft_ref = self.prepare_ref_seq(
                fn=str(input_folder / "hsa_GRCh38_chr8_127735500-127736500.fasta.gz"),
                species="Homo sapiens",
                build="GRCh38",
            )

            host_ref = self.prepare_ref_seq(
                fn=str(input_folder / "mmu_GRCm38_chr15_61857500-61858500.fasta.gz"),
                species="Mus musculus",
                build="GRCm38",
            )

            reads = self.prepare_reads(
                fn=[str(input_folder / "SRR9130497_100_1.fastq.gz")],
            )

            paired_reads = self.prepare_paired_reads(
                mate1=[str(input_folder / "SRR9130497_100_1.fastq.gz")],
                mate2=[str(input_folder / "SRR9130497_100_2.fastq.gz")],
            )

        index = self.run_process(
            process_slug="xengsort-index",
            input_={
                "graft_refs": [graft_ref.id],
                "host_refs": [host_ref.id],
            },
        )

        self.assertFields(index, "graft_species", "Homo sapiens")
        self.assertFields(index, "graft_build", "GRCh38")
        self.assertFields(index, "host_species", "Mus musculus")
        self.assertFields(index, "host_build", "GRCm38")

        self.assertFile(
            index,
            "stats",
            str(output_folder / "index_stats.txt"),
            file_filter=filter_variable_lines,
        )

        self.assertDir(index, "index", output_folder / "xgensort_index.tar.gz")

        classify = self.run_process(
            process_slug="xengsort-classify",
            input_={
                "reads": reads.id,
                "index": index.id,
                "upload_reads": "graft, host",
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFields(classify, "graft_species", "Homo sapiens")
        self.assertFields(classify, "graft_build", "GRCh38")
        self.assertFields(classify, "host_species", "Mus musculus")
        self.assertFields(classify, "host_build", "GRCm38")

        self.assertFile(
            classify,
            "stats",
            str(output_folder / "classification_stats_single.txt"),
            file_filter=filter_variable_lines,
        )

        self.assertFile(
            classify,
            "graft1",
            str(output_folder / "graft_reads_single.fastq.gz"),
            compression="gzip",
        )

        fastq = Data.objects.get(
            process__slug="upload-fastq-single", name="SRR9130497_100_1-graft.fastq.gz"
        )

        self.assertFiles(
            fastq,
            "fastq",
            [str(output_folder / "graft_reads_single.fastq.gz")],
            compression="gzip",
        )

        classify_paired = self.run_process(
            process_slug="xengsort-classify",
            input_={
                "reads": paired_reads.id,
                "index": index.id,
                "upload_reads": "graft",
                "merge_both": True,
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        self.assertFields(classify_paired, "graft_species", "Homo sapiens")
        self.assertFields(classify_paired, "graft_build", "GRCh38")
        self.assertFields(classify_paired, "host_species", "Mus musculus")
        self.assertFields(classify_paired, "host_build", "GRCm38")

        self.assertFile(
            classify_paired,
            "stats",
            str(output_folder / "classification_stats_paired.txt"),
            file_filter=filter_variable_lines,
        )

        self.assertFile(
            classify_paired,
            "graft1",
            str(output_folder / "graft_reads_1.fastq.gz"),
            compression="gzip",
        )

        self.assertFile(
            classify_paired,
            "graft2",
            str(output_folder / "graft_reads_2.fastq.gz"),
            compression="gzip",
        )

        fastq_paired = Data.objects.get(
            process__slug="upload-fastq-paired",
            name="SRR9130497_100_1-graft-both.1.fastq.gz",
        )

        self.assertFiles(
            fastq_paired,
            "fastq",
            [str(output_folder / "graft_reads_1.fastq.gz")],
            compression="gzip",
        )
        self.assertFiles(
            fastq_paired,
            "fastq2",
            [str(output_folder / "graft_reads_2.fastq.gz")],
            compression="gzip",
        )

    @tag_process("rnaseq-vc-preprocess")
    def test_rnaseq_vc_preprocess(self):
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
            star = self.run_process(
                "alignment-star",
                {
                    "reads": reads.id,
                    "genome": star_index.id,
                    "two_pass_mapping": {"two_pass_mode": True},
                    "output_options": {"out_unmapped": True},
                },
            )

        def filter_startedon(line):
            return line.startswith(b"# Started on:") or line.startswith(
                b"# MarkDuplicates"
            )

        inputs = {
            "bam": star.id,
            "ref_seq": ref_seq.id,
            "known_sites": [dbsnp.id],
        }

        preprocess = self.run_process("rnaseq-vc-preprocess", inputs)

        self.assertFile(
            preprocess, "stats", output_folder / "chr1_19000_R1.bam_stats.txt"
        )
        self.assertFile(
            preprocess,
            "metrics_file",
            output_folder / "chr1_19000_R1_markduplicates_metrics.txt",
            file_filter=filter_startedon,
        )
