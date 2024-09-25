from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class GeneFusionWorkflowTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("gene-fusion-calling-arriba")
    def test_gene_fusion_workflow_with_bam(self):
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

        workflow_inputs = {
            "bam": bam.id,
            "gtf": gtf_file.id,
            "genome": genome_file.id,
            "known_fusions_file": known_fusions_file.id,
        }

        self.run_process("gene-fusion-calling-arriba", workflow_inputs)
        fusion_result = Data.objects.filter(process__slug="arriba")[0]

        self.assertFile(
            fusion_result,
            "fusions",
            output_folder / "expected_fusions.tsv.gz",
            compression="gzip",
        )
        self.assertFile(
            fusion_result,
            "discarded_fusions",
            output_folder / "expected_discarded_fusions.tsv.gz",
            compression="gzip",
        )

    @with_resolwe_host
    @tag_process("gene-fusion-calling-arriba")
    def test_gene_fusion_workflow_with_fastq(self):
        input_folder = Path("arriba") / "input"
        output_folder = Path("arriba") / "output"

        with self.preparation_stage():

            reads = self.prepare_paired_reads(
                [input_folder / "reads_R1.fastq.gz"],
                [input_folder / "reads_R2.fastq.gz"],
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

            star_index = self.run_process(
                process_slug="alignment-star-index",
                input_={"ref_seq": genome_file.id, "annotation": gtf_file.id},
            )

            known_fusions_file = self.run_process(
                process_slug="upload-file",
                input_={
                    "src": input_folder / "known_fusions.tsv.gz",
                },
            )

        workflow_inputs = {
            "reads": reads.id,
            "gtf": gtf_file.id,
            "genome": genome_file.id,
            "star_index": star_index.id,
            "known_fusions_file": known_fusions_file.id,
        }

        self.run_process("gene-fusion-calling-arriba", workflow_inputs)
        fusion_result = Data.objects.filter(process__slug="arriba").last()

        self.assertFile(
            fusion_result,
            "fusions",
            output_folder / "expected_fusions.tsv.gz",
            compression="gzip",
        )
        self.assertFile(
            fusion_result,
            "discarded_fusions",
            output_folder / "expected_discarded_fusions.tsv.gz",
            compression="gzip",
        )
