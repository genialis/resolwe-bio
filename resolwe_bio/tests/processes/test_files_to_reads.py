from pathlib import Path

from resolwe.flow.models import Secret
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class FilesToReadsTestCase(BioProcessTestCase):
    @tag_process("basespace-file-import", "files-to-fastq-single")
    def external_test_files_fq_single(self):
        """The following access token was created by the user Jan Otonicar on December
        the 15th, 2021. All the files in this test were uploaded to his basespace profile.
        To create a new acces token follow this link
        https://developer.basespace.illumina.com/docs/content/documentation/authentication/obtaining-access-tokens
        """
        output_folder = Path("basespace_import") / "output"
        with self.preparation_stage():
            # Token with limited scope pre-obtained from dedicated BaseSpace testing app.
            handle = Secret.objects.create_secret(
                "d23ee29d4db4454480ad89d925475913", self.admin
            )

            file_id1 = "25157957505"
            file_id2 = "25157957506"

            import_inputs_1 = {
                "file_id": file_id1,
                "access_token_secret": {"handle": handle},
            }
            basespace_import_1 = self.run_process(
                "basespace-file-import", import_inputs_1
            )

            import_inputs_2 = {
                "file_id": file_id2,
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
            obj=reads,
            field_path="fastq",
            fn_list=[
                output_folder / "single_S1_L001_R1_001.fastq.gz",
                output_folder / "single_S1_L002_R1_001.fastq.gz",
            ],
            compression="gzip",
        )
        del reads.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        del reads.output["fastqc_url"][1]["total_size"]  # Non-deterministic output.
        self.assertFields(
            obj=reads,
            path="fastqc_url",
            value=[
                {
                    "file": "fastqc/test1_S1_L001_R1_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/test1_S1_L001_R1_001_fastqc"],
                },
                {
                    "file": "fastqc/test1_S1_L002_R1_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/test1_S1_L002_R1_001_fastqc"],
                },
            ],
        )

    @tag_process("basespace-file-import", "files-to-fastq-paired")
    def external_test_files_fq_paired(self):
        """The following access token was created by the user Jan Otonicar on December
        the 15th, 2021. All the files in this test were uploaded to his basespace profile.
        To create a new acces token follow this link
        https://developer.basespace.illumina.com/docs/content/documentation/authentication/obtaining-access-tokens
        """
        output_folder = Path("basespace_import") / "output"
        with self.preparation_stage():
            # Token with limited scope pre-obtained from dedicated BaseSpace testing app.
            handle = Secret.objects.create_secret(
                "d23ee29d4db4454480ad89d925475913", self.admin
            )

            file_idU1 = "25158080061"
            file_idU2 = "25158080063"
            file_idD1 = "25158080064"
            file_idD2 = "25158080062"

            upstream_inputs_1 = {
                "file_id": file_idU1,
                "access_token_secret": {"handle": handle},
            }
            basespace_import_upstream_1 = self.run_process(
                "basespace-file-import", upstream_inputs_1
            )

            upstream_inputs_2 = {
                "file_id": file_idU2,
                "access_token_secret": {"handle": handle},
            }
            basespace_import_upstream_2 = self.run_process(
                "basespace-file-import", upstream_inputs_2
            )

            downstream_inputs_1 = {
                "file_id": file_idD1,
                "access_token_secret": {"handle": handle},
            }
            basespace_import_downstream_1 = self.run_process(
                "basespace-file-import", downstream_inputs_1
            )

            downstream_inputs_2 = {
                "file_id": file_idD2,
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
            obj=reads,
            field_path="fastq",
            fn_list=[
                output_folder / "paired_S1_L001_R1_001.fastq.gz",
                output_folder / "paired_S1_L002_R1_001.fastq.gz",
            ],
            compression="gzip",
        )
        self.assertFiles(
            obj=reads,
            field_path="fastq2",
            fn_list=[
                output_folder / "paired_S1_L001_R2_001.fastq.gz",
                output_folder / "paired_S1_L002_R2_001.fastq.gz",
            ],
            compression="gzip",
        )
        del reads.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        del reads.output["fastqc_url"][1]["total_size"]  # Non-deterministic output.
        self.assertFields(
            obj=reads,
            path="fastqc_url",
            value=[
                {
                    "file": "fastqc/test2_S1_L001_R1_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/test2_S1_L001_R1_001_fastqc"],
                },
                {
                    "file": "fastqc/test2_S1_L002_R1_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/test2_S1_L002_R1_001_fastqc"],
                },
            ],
        )
        del reads.output["fastqc_url2"][0]["total_size"]  # Non-deterministic output.
        del reads.output["fastqc_url2"][1]["total_size"]  # Non-deterministic output.
        self.assertFields(
            obj=reads,
            path="fastqc_url2",
            value=[
                {
                    "file": "fastqc/test2_S1_L001_R2_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/test2_S1_L001_R2_001_fastqc"],
                },
                {
                    "file": "fastqc/test2_S1_L002_R2_001_fastqc/fastqc_report.html",
                    "refs": ["fastqc/test2_S1_L002_R2_001_fastqc"],
                },
            ],
        )

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
