from resolwe.flow.models import Secret
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class FilesToReadsTestCase(BioProcessTestCase):
    @tag_process("basespace-file-import", "files-to-fastq-single")
    def external_test_files_fq_single(self):
        with self.preparation_stage():
            # Token with limited scope pre-obtained from dedicated BaseSpace testing app.
            handle = Secret.objects.create_secret(
                "9bdf059c759a429f8af52ca084130060", self.admin
            )

            import_inputs_1 = {
                "file_id": "9461130722",
                "access_token_secret": {"handle": handle},
            }
            basespace_import_1 = self.run_process(
                "basespace-file-import", import_inputs_1
            )

            import_inputs_2 = {
                "file_id": "9461121664",
                "access_token_secret": {"handle": handle},
            }
            basespace_import_2 = self.run_process(
                "basespace-file-import", import_inputs_2
            )

        reads = self.run_process(
            "files-to-fastq-single",
            {"src": [basespace_import_1.pk, basespace_import_2.pk]},
        )

        self.assertFiles(
            reads,
            "fastq",
            ["Test_S1_L001_R1_001.fastq.gz", "Test_S1_L002_R1_001.fastq.gz"],
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

    @tag_process("basespace-file-import", "files-to-fastq-paired")
    def external_test_files_fq_paired(self):
        with self.preparation_stage():
            # Token with limited scope pre-obtained from dedicated BaseSpace testing app.
            handle = Secret.objects.create_secret(
                "d0728b8cceb7455786665453d28c7ebc", self.admin
            )

            upstream_inputs_1 = {
                "file_id": "9864012319",
                "access_token_secret": {"handle": handle},
            }
            basespace_import_upstream_1 = self.run_process(
                "basespace-file-import", upstream_inputs_1
            )

            upstream_inputs_2 = {
                "file_id": "9863993106",
                "access_token_secret": {"handle": handle},
            }
            basespace_import_upstream_2 = self.run_process(
                "basespace-file-import", upstream_inputs_2
            )

            downstream_inputs_1 = {
                "file_id": "9863999826",
                "access_token_secret": {"handle": handle},
            }
            basespace_import_downstream_1 = self.run_process(
                "basespace-file-import", downstream_inputs_1
            )

            downstream_inputs_2 = {
                "file_id": "9863993107",
                "access_token_secret": {"handle": handle},
            }
            basespace_import_downstream_2 = self.run_process(
                "basespace-file-import", downstream_inputs_2
            )

        reads = self.run_process(
            "files-to-fastq-paired",
            {
                "src1": [
                    basespace_import_upstream_1.pk,
                    basespace_import_upstream_2.pk,
                ],
                "src2": [
                    basespace_import_downstream_1.pk,
                    basespace_import_downstream_2.pk,
                ],
            },
        )

        self.assertFiles(
            reads,
            "fastq",
            [
                "TestPaired_S1_L001_R1_001.fastq.gz",
                "TestPaired_S1_L002_R1_001.fastq.gz",
            ],
            compression="gzip",
        )
        self.assertFiles(
            reads,
            "fastq2",
            [
                "TestPaired_S1_L001_R2_001.fastq.gz",
                "TestPaired_S1_L002_R2_001.fastq.gz",
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

    @tag_process("upload-file", "files-to-fastq-single")
    def test_files_fq_single(self):
        with self.preparation_stage():
            lane_1 = self.run_process(
                "upload-file", {"src": "Test_S1_L001_R1_001.fastq.gz"}
            )
            lane_2 = self.run_process(
                "upload-file", {"src": "Test_S1_L002_R1_001.fastq.gz"}
            )

        reads = self.run_process(
            "files-to-fastq-single", {"src": [lane_1.id, lane_2.id]}
        )

        self.assertFiles(
            reads,
            "fastq",
            ["Test_S1_L001_R1_001.fastq.gz", "Test_S1_L002_R1_001.fastq.gz"],
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
            reads, "fastq", ["Test_S1_L001_R1_001_merged.fastq.gz"], compression="gzip"
        )

    @tag_process("upload-file", "files-to-fastq-paired")
    def test_files_fq_paired(self):
        with self.preparation_stage():
            mate_1_lane_1 = self.run_process(
                "upload-file", {"src": "TestPaired_S1_L001_R1_001.fastq.gz"}
            )
            mate_1_lane_2 = self.run_process(
                "upload-file", {"src": "TestPaired_S1_L002_R1_001.fastq.gz"}
            )
            mate_2_lane_1 = self.run_process(
                "upload-file", {"src": "TestPaired_S1_L001_R2_001.fastq.gz"}
            )
            mate_2_lane_2 = self.run_process(
                "upload-file", {"src": "TestPaired_S1_L002_R2_001.fastq.gz"}
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
                "TestPaired_S1_L001_R1_001.fastq.gz",
                "TestPaired_S1_L002_R1_001.fastq.gz",
            ],
            compression="gzip",
        )
        self.assertFiles(
            reads,
            "fastq2",
            [
                "TestPaired_S1_L001_R2_001.fastq.gz",
                "TestPaired_S1_L002_R2_001.fastq.gz",
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
            ["TestPaired_S1_L001_R1_001_merged.fastq.gz"],
            compression="gzip",
        )
        self.assertFiles(
            reads,
            "fastq2",
            ["TestPaired_S1_L001_R2_001_merged.fastq.gz"],
            compression="gzip",
        )
