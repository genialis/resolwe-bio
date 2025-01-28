from pathlib import Path

from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class FilesToReadsTestCase(BioProcessTestCase):
    @tag_process("upload-file", "files-to-fastq-single")
    def test_files_fq_single(self):
        input_folder = Path("test_fastq_upload") / "input"
        output_folder = Path("test_fastq_upload") / "output"
        with self.preparation_stage():
            lane_1 = self.run_process(
                "upload-file", {"src": input_folder / "Test_S1_L001_R1_001.fastq.gz"}
            )
            lane_2 = self.run_process(
                "upload-file", {"src": input_folder / "Test_S1_L002_R1_001.fastq.gz"}
            )

        reads = self.run_process(
            "files-to-fastq-single", {"src": [lane_1.id, lane_2.id]}
        )

        self.assertFiles(
            reads,
            "fastq",
            [
                input_folder / "Test_S1_L001_R1_001.fastq.gz",
                input_folder / "Test_S1_L002_R1_001.fastq.gz",
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
                    "file": "fastqc/Test_S1_L001_R1_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/Test_S1_L001_R1_001_fastqc"],
                },
                {
                    "file": "fastqc/Test_S1_L002_R1_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/Test_S1_L002_R1_001_fastqc"],
                },
            ],
        )

        reads = self.run_process(
            "files-to-fastq-single",
            {
                "src": [lane_1.id, lane_2.id],
                "merge_lanes": True,
            },
        )
        self.assertFiles(
            reads,
            "fastq",
            [output_folder / "Test_S1_L001_R1_001_merged.fastq.gz"],
            compression="gzip",
        )

    @tag_process("upload-file", "files-to-fastq-paired")
    def test_files_fq_paired(self):
        input_folder = Path("test_fastq_upload") / "input"
        output_folder = Path("test_fastq_upload") / "output"
        with self.preparation_stage():
            mate_1_lane_1 = self.run_process(
                "upload-file",
                {"src": input_folder / "TestPaired_S1_L001_R1_001.fastq.gz"},
            )
            mate_1_lane_2 = self.run_process(
                "upload-file",
                {"src": input_folder / "TestPaired_S1_L002_R1_001.fastq.gz"},
            )
            mate_2_lane_1 = self.run_process(
                "upload-file",
                {"src": input_folder / "TestPaired_S1_L001_R2_001.fastq.gz"},
            )
            mate_2_lane_2 = self.run_process(
                "upload-file",
                {"src": input_folder / "TestPaired_S1_L002_R2_001.fastq.gz"},
            )

        reads = self.run_process(
            "files-to-fastq-paired",
            {
                "src1": [mate_1_lane_1.id, mate_1_lane_2.pk],
                "src2": [mate_2_lane_1.pk, mate_2_lane_2.pk],
            },
        )

        self.assertFiles(
            reads,
            "fastq",
            [
                input_folder / "TestPaired_S1_L001_R1_001.fastq.gz",
                input_folder / "TestPaired_S1_L002_R1_001.fastq.gz",
            ],
            compression="gzip",
        )
        self.assertFiles(
            reads,
            "fastq2",
            [
                input_folder / "TestPaired_S1_L001_R2_001.fastq.gz",
                input_folder / "TestPaired_S1_L002_R2_001.fastq.gz",
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
                    "file": "fastqc/TestPaired_S1_L001_R1_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/TestPaired_S1_L001_R1_001_fastqc"],
                },
                {
                    "file": "fastqc/TestPaired_S1_L002_R1_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/TestPaired_S1_L002_R1_001_fastqc"],
                },
            ],
        )
        del reads.output["fastqc_url2"][0]["total_size"]  # Non-deterministic output.
        del reads.output["fastqc_url2"][1]["total_size"]  # Non-deterministic output.
        self.assertFields(
            reads,
            "fastqc_url2",
            [
                {
                    "file": "fastqc/TestPaired_S1_L001_R2_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/TestPaired_S1_L001_R2_001_fastqc"],
                },
                {
                    "file": "fastqc/TestPaired_S1_L002_R2_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/TestPaired_S1_L002_R2_001_fastqc"],
                },
            ],
        )

        reads = self.run_process(
            "files-to-fastq-paired",
            {
                "src1": [mate_1_lane_1.id, mate_1_lane_2.pk],
                "src2": [mate_2_lane_1.pk, mate_2_lane_2.pk],
                "merge_lanes": True,
            },
        )
        self.assertFiles(
            reads,
            "fastq",
            [output_folder / "TestPaired_S1_L001_R1_001_merged.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            reads,
            "fastq2",
            [output_folder / "TestPaired_S1_L001_R2_001_merged.fastq.gz"],
            compression="gzip",
        )
