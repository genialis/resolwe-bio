from pathlib import Path

from django.test import LiveServerTestCase

from resolwe.flow.models import Collection, Data
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.utils.test import BioProcessTestCase


class GeoImportTestCase(BioProcessTestCase, LiveServerTestCase):
    @with_resolwe_host
    @tag_process("geo-import")
    def test_dss_geo(self):
        base = Path("geo_import")
        outputs = base / "outputs"

        collection = Collection.objects.create(
            name="Test collection", contributor=self.contributor
        )

        dss_inputs = {
            "gse_accession": "GSE166144",
            "advanced": {"max_spot_id": 1, "prefetch": False},
        }

        self.run_process("geo-import", dss_inputs, collection=collection)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        sra = Data.objects.filter(process__slug="import-sra-single").last()
        self.assertFiles(
            sra,
            "fastq",
            [str(outputs / "SRR13627912.fastq.gz")],
            compression="gzip",
        )
        metadata = Data.objects.filter(process__slug="upload-metadata-unique").last()
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
            "high throughput sequencing, Expression profiling by array, Genome binding/occupancy profiling by "
            "high throughput sequencing are supported."
        ]
        self.assertEqual(rtpcr.process_error, type_error)

    @with_resolwe_host
    @tag_process("geo-import")
    def test_geo_microarray(self):
        base = Path("geo_import")
        inputs = base / "inputs"
        outputs = base / "outputs"

        collection = Collection.objects.create(
            name="Test collection", contributor=self.contributor
        )

        # Due to Biomart instability the testing is done with the mapping file.
        # For automatic mapping one can comment the mapping_file line.
        dss_inputs = {
            "gse_accession": "GSE172293",
            "advanced": {
                "mapping_file": str(inputs / "HuGene-2_0-st_mapping.tsv.gz"),
                "source": "ENSEMBL",
                "build": "GRCh38.p13",
            },
        }

        self.run_process("geo-import", dss_inputs, collection=collection)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        exp = Data.objects.filter(process__slug="mapped-microarray-expression").last()
        self.assertFile(
            exp,
            "exp",
            str(outputs / "GSM5252014_mapped_exp.tsv.gz"),
            compression="gzip",
        )
        self.assertFields(exp, "species", "Homo sapiens")
        self.assertFields(exp, "source", "ENSEMBL")
        self.assertFields(exp, "build", "GRCh38.p13")
        self.assertFields(
            exp,
            "platform",
            "[HuGene-2_0-st] Affymetrix Human Gene 2.0 ST Array [transcript (gene) version]",
        )
        self.assertFields(exp, "platform_id", "GPL16686")
        self.assertFields(exp, "probe_mapping", "Custom")

        metadata = Data.objects.filter(process__slug="upload-metadata-unique").last()
        self.assertFile(
            metadata,
            "table",
            str(outputs / "GSE172293_metadata.tsv"),
        )

    @with_resolwe_host
    @tag_process("geo-import")
    def test_geo_chipseq(self):
        base = Path("geo_import")
        outputs = base / "outputs"

        collection = Collection.objects.create(
            name="Test collection", contributor=self.contributor
        )

        dss_inputs = {
            "gse_accession": "GSE176232",
            "advanced": {"max_spot_id": 1, "prefetch": False},
        }

        self.run_process("geo-import", dss_inputs, collection=collection)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        sra = Data.objects.filter(process__slug="import-sra-single").last()
        self.assertFiles(
            sra,
            "fastq",
            [str(outputs / "SRR14743655.fastq.gz")],
            compression="gzip",
        )
        metadata = Data.objects.filter(process__slug="upload-metadata-unique").last()
        self.assertFile(
            metadata,
            "table",
            str(outputs / "GSE176232_metadata.tsv"),
        )
        del sra.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        self.assertFields(
            sra,
            "fastqc_url",
            [
                {
                    "file": "fastqc/SRR14743655_fastqc/fastqc_report.html",
                    "refs": ["fastqc/SRR14743655_fastqc"],
                },
            ],
        )
        descriptor = {
            "general": {
                "species": "Homo sapiens",
                "annotator": "Matthew,T,Weirauch",
                "cell_line": "GM12878 (B-Lymphocyte) LCL",
                "description": "ChIP_EBNA2_GM12878_rep2-E00457",
                "biosample_type": "cell_line",
                "growth_protocol": "Cells were cultured in 10% FBS supplemented RPMI 1640 medium for 2 weeks.",
                "biosample_source": "GM12878 (B-Lymphocyte) LCL",
                "treatment_protocol": "",
            },
            "experiment": {
                "molecule": "genomic_dna",
                "assay_type": "chip-seq",
                "extract_protocol": (
                    "Cells were crosslinked and nuclei were sonicated as described previously (Lu et al. 2015). "
                    "Libraries were prepared via ChIPmentation (Schmidl et al. 2015)."
                ),
            },
        }

        self.assertEqual(sra.entity.descriptor, descriptor)

    @with_resolwe_host
    @tag_process("geo-import")
    def test_geo_ena(self):
        base = Path("geo_import")
        outputs = base / "outputs"

        collection = Collection.objects.create(
            name="Test collection", contributor=self.contributor
        )

        dss_inputs = {
            "gse_accession": "GSE24998",
            "advanced": {"max_spot_id": 1, "prefetch": False},
        }

        self.run_process("geo-import", dss_inputs, collection=collection)

        for data in Data.objects.all():
            self.assertStatus(data, Data.STATUS_DONE)

        metadata = Data.objects.filter(process__slug="upload-metadata-unique").last()
        self.assertFile(
            metadata,
            "table",
            str(outputs / "GSE24998_metadata.tsv"),
        )

        sra = Data.objects.filter(process__slug="import-sra-single").last()
        self.assertFiles(
            sra,
            "fastq",
            [str(outputs / "ERR019653.fastq.gz")],
            compression="gzip",
        )

        del sra.output["fastqc_url"][0]["total_size"]  # Non-deterministic output.
        self.assertFields(
            sra,
            "fastqc_url",
            [
                {
                    "file": "fastqc/ERR019653_fastqc/fastqc_report.html",
                    "refs": ["fastqc/ERR019653_fastqc"],
                },
            ],
        )

        sample_descriptor = {
            "general": {
                "annotator": "ArrayExpress EBI",
                "description": "Provider: EMBL_Heidelberg",
                "growth_protocol": "grow | The E. coli K-12 MG1655 bacterial strains used in this "
                "work are the following: E. coli MG1655 (F- lambda- ilvG- rfb-50 rph-1). "
                "Luria-Bertani (0.5% NaCl) broth and agar (15 g/liter) were used for routine "
                "growth.",
                "biosample_source": "E_coli_K12_MG1655_RNAseq_LB_Transition_to_stationary",
            },
            "experiment": {
                "molecule": "total_rna",
                "assay_type": "rna-seq",
                "extract_protocol": "nucleic_acid_extraction | To prepare cells for RNA "
                "extraction, 100 ml of fresh LB was inoculated 1:200 from an overnight culture "
                "in a 250 ml flask and incubated with shaking at 180 r.p.m. in a New Brunswick "
                "C76 waterbath at 37C. Two biological replicates were performed for each strain "
                "and samples were taken at early-exponential, mid-exponential, "
                "transition-to-stationary and stationary phase. The cells were pelleted by "
                "centrifugation (10000 g, 10 min, 4¡C), washed in 1xPBS and pellets were "
                "snap-frozen and stored at -80C until required. RNA was extracted using Trizol "
                "Reagent (Invitrogen) according to the manufacturer's protocol until the "
                "chloroform extraction step. The aqueous phase was then loaded onto mirVanaTM "
                "miRNA Isolation kit (Ambion Inc.) columns and washed according to the "
                "manufacturer's protocol. Total RNA was eluted in 50µl of RNAase free water. "
                "The concentration was then determined using a Nanodrop ND-1000 machine (NanoDrop "
                "Technologies), and RNA quality was tested by visualization on agarose gels and "
                "by Agilent 2100 Bioanalyser (Agilent Technologies). sequencing | Standard "
                "Illumina protocol for cDNA sequencing",
            },
        }

        self.assertEqual(sra.entity.descriptor, sample_descriptor)
