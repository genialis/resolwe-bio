from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase


class SamtoolsProcessorTestCase(BioProcessTestCase):
    @tag_process("samtools-view")
    def samtools_view(self):
        input_folder = Path("samtools") / "inputs"
        output_folder = Path("samtools") / "outputs"
        with self.preparation_stage():
            inputs_bam = {
                "src": input_folder / "samtools_in.bam",
                "species": "Homo sapiens",
                "build": "GRCh38",
            }
            bam = self.run_process("upload-bam", inputs_bam)

            inputs_bed = {
                "src": input_folder / "regions_bed.bed",
                "species": "Homo sapiens",
                "build": "GRCh38",
            }
            bed = self.run_process("upload-bed", inputs_bed)

        inputs = {
            "bam": bam.id,
            "region": "7:116755355-116755480",
        }
        samtools = self.run_process("samtools-view", inputs)
        self.assertFile(samtools, "stats", output_folder / "regions.bam_stats.txt")

        inputs = {
            "bam": bam.id,
            "bedfile": bed.id,
            "advanced": {
                "include_header": False,
            },
        }
        samtools = self.run_process("samtools-view", inputs)
        self.assertFile(samtools, "stats", output_folder / "bedfile.bam_stats.txt")
        self.assertFile(samtools, "bam", output_folder / "samtools_bed.bam")

        inputs = {
            "bam": bam.id,
            "bedfile": bed.id,
            "advanced": {"subsample": 0.8, "subsample_seed": 60},
        }
        samtools = self.run_process("samtools-view", inputs)
        self.assertFileExists(samtools, "stats")

        inputs = {
            "bam": bam.id,
            "region": "2:50493-50634",
            "advanced": {"only_header": True},
        }
        samtools = self.run_process("samtools-view", inputs)
        self.assertFileExists(samtools, "bam")

        inputs = {
            "bam": bam.id,
        }
        samtools = self.run_process("samtools-view", inputs, Data.STATUS_ERROR)
        error_msg = [("No region or BED file specified.")]
        self.assertEqual(samtools.process_error, error_msg)

        inputs = {
            "bam": bam.id,
            "region": "chr7:116755355-116755480",
        }
        samtools = self.run_process("samtools-view", inputs)
        warning_msg = [
            (
                '[main_samview] region "chr7:116755355-116755480" '
                "specifies an invalid region or unknown reference. Continue anyway.\n"
            )
        ]
        self.assertEqual(samtools.process_warning, warning_msg)

    @tag_process("samtools-coverage-single", "samtools-coverage-multi")
    def samtools_coverage(self):
        input_folder = Path("samtools") / "inputs"
        output_folder = Path("samtools") / "outputs"
        with self.preparation_stage():
            inputs_bam = {
                "src": input_folder / "samtools_in.bam",
                "species": "Homo sapiens",
                "build": "GRCh38",
            }
            bam1 = self.run_process("upload-bam", inputs_bam)

            inputs_bam = {
                "src": output_folder / "samtools_bed.bam",
                "species": "Homo sapiens",
                "build": "GRCh38",
            }
            bam2 = self.run_process("upload-bam", inputs_bam)

            inputs_bam = {
                "src": output_folder / "samtools_bed.bam",
                "species": "Homo sapiens",
                "build": "GRCh37",
            }
            bam3 = self.run_process("upload-bam", inputs_bam)

        inputs = {
            "bam": [bam1.id, bam2.id],
            "region": "7",
            "advanced": {"min_read_length": 27, "depth": 0},
        }
        coverage = self.run_process("samtools-coverage-multi", inputs)
        self.assertFile(coverage, "table", output_folder / "coverage_multi.tsv")

        inputs = {
            "bam": [bam1.id, bam2.id, bam3.id],
            "region": "7",
            "advanced": {"min_read_length": 27, "depth": 0},
        }
        coverage = self.run_process(
            "samtools-coverage-multi", inputs, Data.STATUS_ERROR
        )
        error_msg = [
            (
                "Not all BAM files have the same genome build. "
                "BAM file samtools_in.bam has build GRCh38, while file "
                "samtools_bed.bam has build GRCh37."
            )
        ]
        self.assertEqual(coverage.process_error, error_msg)

        inputs = {
            "bam": bam1.id,
            "region": "7:116755355-116755480",
            "advanced": {
                "excl_flags": ["SECONDARY", "QCFAIL"],
            },
        }
        coverage = self.run_process("samtools-coverage-single", inputs)
        self.assertFile(coverage, "table", output_folder / "coverage_single.tsv")

        inputs = {
            "bam": bam1.id,
            "region": "chr7:116755355-1167554804",
            "advanced": {
                "excl_flags": ["SECONDARY"],
            },
        }
        coverage = self.run_process(
            "samtools-coverage-single", inputs, Data.STATUS_ERROR
        )

    @tag_process("samtools-bedcov")
    def samtools_bedcov(self):
        input_folder = Path("samtools") / "inputs"
        output_folder = Path("samtools") / "outputs"
        with self.preparation_stage():
            inputs_bam = {
                "src": input_folder / "samtools_in.bam",
                "species": "Homo sapiens",
                "build": "GRCh38",
            }
            bam = self.run_process("upload-bam", inputs_bam)

            inputs_bed = {
                "src": input_folder / "regions_bed.bed",
                "species": "Homo sapiens",
                "build": "GRCh38",
            }
            bed = self.run_process("upload-bed", inputs_bed)

        inputs = {"bam": bam.id, "bedfile": bed.id}
        samtools = self.run_process("samtools-bedcov", inputs)
        self.assertFile(
            samtools,
            "coverage_report",
            output_folder / "regions_bed_coverage.txt.gz",
            compression="gzip",
        )
        inputs = {
            "bam": bam.id,
            "bedfile": bed.id,
            "advanced": {
                "output_option": "mean",
                "min_read_qual": 5,
            },
        }
        samtools = self.run_process("samtools-bedcov", inputs)
        self.assertFile(
            samtools,
            "coverage_report",
            output_folder / "regions_bed_coverage_mean.txt.gz",
            compression="gzip",
        )

    @tag_process("samtools-depth-single")
    def samtools_depth_single(self):
        input_folder = Path("samtools") / "inputs"
        output_folder = Path("samtools") / "outputs"
        with self.preparation_stage():
            inputs_bam = {
                "src": input_folder / "samtools_in.bam",
                "species": "Homo sapiens",
                "build": "GRCh38",
            }
            bam = self.run_process("upload-bam", inputs_bam)

            inputs_bed = {
                "src": input_folder / "regions_bed.bed",
                "species": "Homo sapiens",
                "build": "GRCh38",
            }
            bed = self.run_process("upload-bed", inputs_bed)

        inputs = {"bam": bam.id, "bedfile": bed.id}
        samtools = self.run_process("samtools-depth-single", inputs)
        self.assertFile(
            obj=samtools,
            field_path="depth_report",
            fn=output_folder / "samtools_in_depth.txt.gz",
            compression="gzip",
        )
