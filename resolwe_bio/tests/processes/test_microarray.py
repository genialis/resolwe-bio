from pathlib import Path

from django.core.management import call_command

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import (
    TEST_FILES_DIR,
    KBBioProcessTestCase,
    skipUnlessLargeFiles,
)


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


class MethylationArrayTestCase(KBBioProcessTestCase):
    @skipUnlessLargeFiles("methylation_beta_values_annotated.txt.gz")
    @tag_process("methylation-array-sesame")
    def test_methylation_array(self):
        base = Path("methylation")
        test_dir = Path(TEST_FILES_DIR)
        inputs = base / "inputs"
        outputs = base / "output"
        large = Path("large")

        with self.preparation_stage():
            probe_mappings = test_dir / inputs / "illumina_probe_ids_to_ensembl.tab"
            call_command("insert_mappings", str(probe_mappings))

            species = "Homo sapiens"
            platform = "HM450"

            idat = self.run_process(
                process_slug="upload-idat",
                input_={
                    "red_channel": large / "6042316072_R03C01_Red.idat.gz",
                    "green_channel": large / "6042316072_R03C01_Grn.idat.gz",
                    "species": species,
                    "platform": platform,
                },
            )

        metarr = self.run_process(
            process_slug="methylation-array-sesame",
            input_={
                "idat_file": idat.id,
            },
        )

        self.assertFields(obj=metarr, path="species", value=species)
        self.assertFields(obj=metarr, path="platform", value=platform)
        self.assertFile(
            obj=metarr,
            field_path="qc_data",
            fn=str(outputs / "QC_data.txt"),
        )

        self.assertFile(
            obj=metarr,
            field_path="methylation_data",
            fn=str(large / "methylation_beta_values_annotated.txt.gz"),
            compression="gzip",
        )
