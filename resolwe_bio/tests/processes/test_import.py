import os

from resolwe.flow.models import Data, Secret
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class ImportProcessorTestCase(BioProcessTestCase):
    @tag_process("import-sra", "import-sra-single", "import-sra-paired")
    def external_test_sra(self):
        # single-end reads from Polyak RNA-seq demo dataset
        inputs = {
            "sra_accession": ["SRR1661332", "SRR1661333"],
            "advanced": {"max_spot_id": 1,},
        }
        self.run_process("import-sra", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        import_sra = Data.objects.last()
        self.assertFiles(
            import_sra,
            "fastq",
            [
                os.path.join("import-sra", "output", "SRR1661332.fastq.gz"),
                os.path.join("import-sra", "output", "SRR1661333.fastq.gz"),
            ],
            compression="gzip",
        )
        del import_sra.output["fastqc_url"][0][
            "total_size"
        ]  # Non-deterministic output.
        del import_sra.output["fastqc_url"][1][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFields(
            import_sra,
            "fastqc_url",
            [
                {
                    "file": "fastqc/SRR1661332_fastqc/fastqc_report.html",
                    "refs": ["fastqc/SRR1661332_fastqc"],
                },
                {
                    "file": "fastqc/SRR1661333_fastqc/fastqc_report.html",
                    "refs": ["fastqc/SRR1661333_fastqc"],
                },
            ],
        )

        # paired-end reads from Zoghbi RNA-seq demo dataset
        inputs = {
            "sra_accession": ["SRR2124780", "SRR2124781"],
            "advanced": {"max_spot_id": 1,},
        }
        self.run_process("import-sra", inputs)
        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)
        import_sra = Data.objects.last()
        self.assertFiles(
            import_sra,
            "fastq",
            [
                os.path.join("import-sra", "output", "SRR2124780.fastq.gz"),
                os.path.join("import-sra", "output", "SRR2124781.fastq.gz"),
            ],
            compression="gzip",
        )
        del import_sra.output["fastqc_url"][0][
            "total_size"
        ]  # Non-deterministic output.
        del import_sra.output["fastqc_url"][1][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFields(
            import_sra,
            "fastqc_url",
            [
                {
                    "file": "fastqc/SRR2124780_1_fastqc/fastqc_report.html",
                    "refs": ["fastqc/SRR2124780_1_fastqc"],
                },
                {
                    "file": "fastqc/SRR2124781_1_fastqc/fastqc_report.html",
                    "refs": ["fastqc/SRR2124781_1_fastqc"],
                },
            ],
        )
        del import_sra.output["fastqc_url2"][0][
            "total_size"
        ]  # Non-deterministic output.
        del import_sra.output["fastqc_url2"][1][
            "total_size"
        ]  # Non-deterministic output.
        self.assertFields(
            import_sra,
            "fastqc_url2",
            [
                {
                    "file": "fastqc/SRR2124780_2_fastqc/fastqc_report.html",
                    "refs": ["fastqc/SRR2124780_2_fastqc"],
                },
                {
                    "file": "fastqc/SRR2124781_2_fastqc/fastqc_report.html",
                    "refs": ["fastqc/SRR2124781_2_fastqc"],
                },
            ],
        )

        # test error messages
        sra = self.run_process(
            "import-sra-single",
            {"sra_accession": ["non-existing SRR"]},
            Data.STATUS_ERROR,
        )
        self.assertEqual(
            sra.process_error[0],
            "Download of non-existing SRR reads with fastq-dump failed.",
        )
        sra = self.run_process(
            "import-sra-paired",
            {"sra_accession": ["non-existing SRR"]},
            Data.STATUS_ERROR,
        )
        self.assertEqual(
            sra.process_error[0],
            "Download of non-existing SRR reads with fastq-dump failed.",
        )

        # run on SRA samples of different types (single-end and paired-end)
        sra = self.run_process(
            "import-sra",
            {"sra_accession": ["SRR1661332", "SRR2124781"]},
            Data.STATUS_ERROR,
        )
        self.assertEqual(
            sra.process_error[0],
            (
                "All reads must be either single-end or paired-end. Mixing SRA samples of "
                "different types is not allowed."
            ),
        )

    @tag_process("basespace-file-import")
    def external_test_basespace_import(self):
        # Token with limited scope pre-obtained from dedicated BaseSpace testing app.
        handle = Secret.objects.create_secret(
            "9bdf059c759a429f8af52ca084130060", self.admin
        )

        inputs = {"file_id": "9461130722", "access_token_secret": {"handle": handle}}
        file = self.run_process("basespace-file-import", inputs)

        self.assertFile(
            file, "file", "Test_S1_L001_R1_001.fastq.gz", compression="gzip"
        )
