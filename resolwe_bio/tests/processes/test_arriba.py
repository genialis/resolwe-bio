from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class ArribaProcessorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("arriba")
    def test_arriba_fusion_detection(self):
        input_folder = Path("arriba") / "input"
        output_folder = Path("arriba") / "output"

        with self.preparation_stage():
            bam = self.run_process(
                process_slug="upload-bam",
                input_={
                    "src": input_folder / "aligned_samples.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            bam_error = self.run_process(
                process_slug="upload-bam",
                input_={
                    "src": input_folder / "aligned_samples_error.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            gtf_file = self.run_process(
                process_slug="upload-gtf",
                input_={
                    "src": input_folder / "minigenome.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            genome_file = self.run_process(
                process_slug="upload-fasta-nucl",
                input_={
                    "src": input_folder / "minigenome.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            known_fusions_file = self.run_process(
                process_slug="upload-file",
                input_={
                    "src": input_folder / "known_fusions.tsv.gz",
                },
            )

        arriba_inputs = {
            "bam": bam.id,
            "gtf": gtf_file.id,
            "genome": genome_file.id,
            "known_fusions_file": known_fusions_file.id,
        }

        arriba = self.run_process("arriba", arriba_inputs)

        self.assertFile(
            arriba,
            "fusions",
            output_folder / "expected_fusions.tsv.gz",
            compression="gzip",
        )
        self.assertFile(
            arriba,
            "discarded_fusions",
            output_folder / "expected_discarded_fusions.tsv.gz",
            compression="gzip",
        )

        arriba_inputs["bam"] = bam_error.id

        arriba = self.run_process("arriba", arriba_inputs, Data.STATUS_ERROR)

        error_msg = [
            (
                "STAR parameters were not set correctly. "
                "Please check that the --chimSegmentMin parameter was passed when aligning using STAR and was not zero."
            )
        ]

        self.assertEqual(arriba.process_error, error_msg)
