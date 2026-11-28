from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class FusionAnnotationProcessorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("fusion-annotation")
    def test_fusion_annotation(self):
        arriba_folder = Path("arriba") / "input"
        input_folder = Path("fusion_annotation") / "input"
        output_folder = Path("fusion_annotation") / "output"

        with self.preparation_stage():
            bam = self.run_process(
                process_slug="upload-bam",
                input_={
                    "src": arriba_folder / "aligned_samples.bam",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            gtf_file = self.run_process(
                process_slug="upload-gtf",
                input_={
                    "src": arriba_folder / "minigenome.gtf.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            genome_file = self.run_process(
                process_slug="upload-fasta-nucl",
                input_={
                    "src": arriba_folder / "minigenome.fasta.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

            known_fusions_file = self.run_process(
                process_slug="upload-file",
                input_={
                    "src": arriba_folder / "known_fusions.tsv.gz",
                },
            )

            protein_domains_file = self.run_process(
                process_slug="upload-file",
                input_={
                    "src": arriba_folder / "protein_domains.gff3.gz",
                },
            )

            arriba = self.run_process(
                "arriba",
                {
                    "bam": bam.id,
                    "gtf": gtf_file.id,
                    "genome": genome_file.id,
                    "known_fusions_file": known_fusions_file.id,
                    "protein_domains_file": protein_domains_file.id,
                },
            )

            snapshot = self.run_process(
                process_slug="upload-file",
                input_={
                    "src": input_folder / "snapshot.json",
                },
            )

        annotation = self.run_process(
            "fusion-annotation",
            {
                "fusions": arriba.id,
                "snapshot": snapshot.id,
            },
        )

        self.assertFile(
            annotation,
            "annotated_fusions",
            output_folder / "expected_annotated_fusions.tsv.gz",
            compression="gzip",
        )
        self.assertFile(
            annotation,
            "summary",
            output_folder / "expected_annotation_summary.json",
        )
        self.assertFields(annotation, "snapshot_version", "test-2026-11-25")

        annotation = self.run_process(
            "fusion-annotation",
            {
                "fusions": arriba.id,
                "snapshot": known_fusions_file.id,
            },
            Data.STATUS_ERROR,
        )

        self.assertEqual(
            annotation.process_error,
            ["The knowledge base snapshot is not a valid JSON file."],
        )
