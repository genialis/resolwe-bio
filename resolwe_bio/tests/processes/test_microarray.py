from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class MicroarrayProcessorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("upload-microarray-expression", "mapped-microarray-expression")
    def test_microarray_upload(self):
        base = Path("upload_microarray")
        inputs = base / "input"
        outputs = base / "output"

        unmapped_path = str(inputs / "GSM971957_unmapped_exp.tsv.gz")

        unmapped = self.run_process(
            "upload-microarray-expression",
            {
                "exp": unmapped_path,
                "exp_type": "log2 normalized intensity signal",
                "platform": "Affymetrix Human Genome U133 Plus 2.0 Array",
                "species": "Homo sapiens",
            },
        )
        self.assertFile(unmapped, "exp", unmapped_path, compression="gzip")
        self.assertFields(unmapped, "exp_type", "log2 normalized intensity signal")
        self.assertFields(
            unmapped, "platform", "Affymetrix Human Genome U133 Plus 2.0 Array"
        )
        self.assertFields(unmapped, "species", "Homo sapiens")

        mapped_path = str(inputs / "GSM971957_mapped_exp.tsv.gz")
        mapped = self.run_process(
            "mapped-microarray-expression",
            {
                "exp_unmapped": unmapped.id,
                "exp": mapped_path,
                "source": "ENSEMBL",
                "build": "GRCh38.p13",
                "probe_mapping": "affy_hg_u133_plus_2",
            },
        )
        self.assertFile(mapped, "exp", mapped_path, compression="gzip")
        self.assertJSON(
            mapped, mapped.output["exp_json"], "", str(outputs / "GSM971958.json.gz")
        )
        self.assertFields(mapped, "exp_type", "log2 normalized intensity signal")
        self.assertFields(
            mapped, "platform", "Affymetrix Human Genome U133 Plus 2.0 Array"
        )
        self.assertFile(
            mapped,
            "exp_set",
            str(outputs / "GSM971958_expressions.txt.gz"),
            compression="gzip",
        )
        self.assertFields(mapped, "source", "ENSEMBL")
        self.assertFields(mapped, "species", "Homo sapiens")
        self.assertFields(mapped, "build", "GRCh38.p13")
        self.assertFields(mapped, "feature_type", "gene")
        self.assertFields(mapped, "probe_mapping", "affy_hg_u133_plus_2")

    @with_resolwe_host
    @tag_process("map-microarray-probes")
    def test_probe_mapping(self):
        base = Path("microarray_mapping")
        inputs = base / "input"
        outputs = base / "output"
        with self.preparation_stage():
            unmapped_path = str(inputs / "GSM971957_unmapped_exp.tsv.gz")
            unmapped = self.run_process(
                "upload-microarray-expression",
                {
                    "exp": unmapped_path,
                    "exp_type": "log2 normalized intensity signal",
                    "platform": "Affymetrix Human Genome U133 Plus 2.0 Array",
                    "platform_id": "GPL570",
                    "species": "Homo sapiens",
                },
            )
            unmapped = self.run_process(
                "upload-microarray-expression",
                {
                    "exp": unmapped_path,
                    "exp_type": "log2 normalized intensity signal",
                    "platform": "Affymetrix Human Genome U133 Plus 2.0 Array",
                    "platform_id": "GPL570",
                    "species": "Homo sapiens",
                },
            )

        # Due to Biomart instability the testing is done with the mapping file.
        # For automatic mapping one can comment the mapping_file line.
        self.run_process(
            "map-microarray-probes",
            {
                "expressions": [unmapped.id],
                "mapping_file": str(inputs / "mapping.tsv"),
                "source": "ENSEMBL",
                "build": "GRCh38.p13",
            },
        )

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        mapped = Data.objects.filter(
            process__slug="mapped-microarray-expression"
        ).last()

        self.assertFile(
            mapped, "exp", str(outputs / "GSM971957_mapped.tsv.gz"), compression="gzip"
        )
        self.assertFields(mapped, "source", "ENSEMBL")
        self.assertFields(mapped, "build", "GRCh38.p13")
        self.assertFields(
            mapped, "platform", "Affymetrix Human Genome U133 Plus 2.0 Array"
        )
        self.assertFields(mapped, "probe_mapping", "Custom")
