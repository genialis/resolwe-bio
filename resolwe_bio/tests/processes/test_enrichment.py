from pathlib import Path

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import KBBioProcessTestCase


class EnrichmentProcessorTestCase(KBBioProcessTestCase):
    @with_resolwe_host
    @tag_process("goenrichment")
    def test_go_enrichment_dicty(self):
        input_folder = Path("goenrichment") / "input"
        output_folder = Path("goenrichment") / "output"
        with self.preparation_stage():
            inputs = {"src": input_folder / "ontology_dicty_cropped.obo.gz"}
            ontology = self.run_process("upload-obo", inputs)

            inputs = {
                "src": input_folder / "gaf_dicty_cropped.gz",
                "source": "DICTYBASE",
                "species": "Dictyostelium discoideum",
            }
            annotation = self.run_process("upload-gaf", inputs)

        inputs = {
            "ontology": ontology.pk,
            "gaf": annotation.pk,
            "pval_threshold": 1,
            "source": "DICTYBASE",
            "species": "Dictyostelium discoideum",
            "genes": ["DDB_G0277589", "DDB_G0286855", "DDB_G0267640"],
        }

        enrichment = self.run_process("goenrichment", inputs)
        self.assertEqual(len(enrichment.process_warning), 0)
        self.assertJSON(
            enrichment,
            enrichment.output["terms"],
            "",
            output_folder / "go_enriched_terms_dicty.json.gz",
        )

        inputs = {
            "ontology": ontology.pk,
            "gaf": annotation.pk,
            "pval_threshold": 1,
            "source": "DICTYBASE",
            "species": "Mus musculus",
            "genes": ["DDB_G0277589", "DDB_G0286855", "DDB_G0267640"],
        }
        enrichment = self.run_process("goenrichment", inputs, Data.STATUS_ERROR)
        warning_msg = [
            "Selected genes Species must be the same as the Species "
            "field of the GAF file."
        ]
        error_msg = [
            "Selected genes are from Mus musculus, while GAF "
            "file has defined Dictyostelium discoideum under Species field."
        ]
        self.assertEqual(enrichment.process_warning, warning_msg)
        self.assertEqual(enrichment.process_error, error_msg)

    @with_resolwe_host
    @tag_process("goenrichment")
    def test_go_enrichment(self):
        input_folder = Path("goenrichment") / "input"
        output_folder = Path("goenrichment") / "output"
        with self.preparation_stage():
            inputs = {"src": input_folder / "go.obo.gz"}
            ontology_hs = self.run_process("upload-obo", inputs)

            inputs = {
                "src": input_folder / "goa_human.gaf.txt.gz",
                "source": "UniProtKB",
                "species": "Homo sapiens",
            }
            gaf_hs = self.run_process("upload-gaf", inputs)
            inputs = {"src": input_folder / "ontology_mus_cropped.obo.gz"}
            ontology = self.run_process("upload-obo", inputs)

            inputs = {
                "src": input_folder / "gaf_mgi_cropped.gz",
                "source": "MGI",
                "species": "Mus musculus",
            }
            gaf = self.run_process("upload-gaf", inputs)

            inputs = {
                "src": input_folder / "gaf_mgi_cropped.gz",
                "source": "AFFY",
                "species": "Mus musculus",
            }
            gaf_error = self.run_process("upload-gaf", inputs)

        inputs = {
            "ontology": ontology_hs.pk,
            "gaf": gaf_hs.pk,
            "genes": [
                "ENSG00000247315",
                "ENSG00000125875",
                "ENSG00000101255",
                "ENSG00000125841",
            ],
            "source": "ENSEMBL",
            "species": "Homo sapiens",
        }
        enrichment = self.run_process("goenrichment", inputs)
        self.assertFile(enrichment, "ids", output_folder / "mapped_ids.txt")
        self.assertEqual(
            enrichment.process_warning[0],
            "Mapping ENSG00000125841 returned multiple times.",
        )
        inputs = {
            "ontology": ontology.pk,
            "gaf": gaf.pk,
            "pval_threshold": 1,
            "genes": ["193202", "56535"],
            "source": "NCBI",
            "species": "Mus musculus",
        }

        enrichment = self.run_process("goenrichment", inputs)
        self.assertEqual(len(enrichment.process_warning), 1)
        self.assertEqual(
            enrichment.process_warning[0], "Not all features could be mapped."
        )
        self.assertJSON(
            enrichment,
            enrichment.output["terms"],
            "",
            output_folder / "go_enriched_terms_mouse.json.gz",
        )

        inputs = {
            "ontology": ontology.pk,
            "gaf": gaf.pk,
            "pval_threshold": 1,
            "genes": ["GCTAG33G", "CGAGCG323"],
            "source": "NCBI",
            "species": "Mus musculus",
        }

        enrichment = self.run_process("goenrichment", inputs, Data.STATUS_ERROR)
        error_msg = ["No genes were fetched from the knowledge base."]
        self.assertEqual(enrichment.process_error, error_msg)

        inputs = {
            "ontology": ontology.pk,
            "gaf": gaf_error.pk,
            "pval_threshold": 1,
            "genes": ["193202", "56535"],
            "source": "NCBI",
            "species": "Mus musculus",
        }

        enrichment = self.run_process("goenrichment", inputs, Data.STATUS_ERROR)
        error_msg = ["Failed to map features."]
        self.assertEqual(enrichment.process_error, error_msg)
