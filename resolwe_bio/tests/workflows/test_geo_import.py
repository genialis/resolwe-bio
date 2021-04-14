from pathlib import Path

from django.test import LiveServerTestCase

from resolwe.flow.models import Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import BioProcessTestCase


class GeoImportTestCase(BioProcessTestCase, LiveServerTestCase):
    @with_resolwe_host
    @tag_process("geo-import")
    def test_dss_geo(self):
        base = Path("geo_import")
        outputs = base / "outputs"
        dss_inputs = {
            "gse_accession": "GSE166144",
            "advanced": {"max_spot_id": 1, "prefetch": False},
        }

        self.run_process("geo-import", dss_inputs)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        sra = Data.objects.filter(process__slug="import-sra-single").last()
        self.assertFiles(
            sra,
            "fastq",
            [str(outputs / "SRR13627912.fastq.gz")],
            compression="gzip",
        )
        metadata = Data.objects.filter(process__slug="upload-orange-metadata").last()
        self.assertFile(
            metadata,
            "table",
            str(outputs / "GSE166144_metadata.tsv"),
        )
        del sra.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        self.assertFields(
            sra,
            "fastqc_url",
            [
                {
                    "file": "fastqc/SRR13627912_fastqc/fastqc_report.html",
                    "refs": ["fastqc/SRR13627912_fastqc"],
                },
            ],
        )
        sra = Data.objects.filter(process__slug="import-sra-single").last()
        descriptor = {
            "general": {
                "species": "Homo sapiens",
                "biosample_type": "cell_line",
                "cell_line": "HepG2",
                "annotator": "Lin He",
                "description": "ING5 knockdown",
                "biosample_source": "HepG2 cells",
                "growth_protocol": "HepG2 cells were maintained according to the ATCC’s recommendation.",
                "treatment_protocol": (
                    "HepG2 cells were transfected with vector or FLAG-JFK or treated with control siRNA or ING5 siRNA."
                ),
            },
            "experiment": {
                "assay_type": "rna-seq",
                "extract_protocol": (
                    "Total mRNAs were isolated with Trizol reagents (Invitrogen) for cDNA synthesis, library "
                    "construction, and sequencing using HiSeq 2500. RNA libraries were prepared for sequencing using "
                    "standard Illumina protocols."
                ),
                "molecule": "total_rna",
                "platform": "hiseq_2000",
            },
        }

        self.assertEqual(sra.entity.descriptor, descriptor)

        # Non-existant GSE series.
        wrong = self.run_process(
            "geo-import", {"gse_accession": "GSE99999999"}, Data.STATUS_ERROR
        )
        gse_error = [
            "Download of GSE99999999 failed. ID could be incorrect or the data might not be public yet."
        ]
        self.assertEqual(wrong.process_error, gse_error)
        # GEO dataset number.
        gds = self.run_process(
            "geo-import", {"gse_accession": "GDS6063"}, Data.STATUS_ERROR
        )
        accsession_error = [
            "GEO series accessions (GSE) are supported but GDS6063 was provided."
        ]
        self.assertEqual(gds.process_error, accsession_error)

        # RT-PCR data.
        rtpcr = self.run_process(
            "geo-import", {"gse_accession": "GSE165150"}, Data.STATUS_ERROR
        )
        type_error = [
            "No supported series types found. Got Expression profiling by RT-PCR but only Expression profiling by "
            "high throughput sequencing are supported."
        ]
        self.assertEqual(rtpcr.process_error, type_error)
