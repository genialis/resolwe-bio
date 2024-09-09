from pathlib import Path

from django.contrib.auth.models import User

from resolwe.flow.models import Data
from resolwe.flow.models import Entity as Sample
from resolwe.test import tag_process

from resolwe_bio.utils.test import BioProcessTestCase
from resolwe_bio.variants.models import Variant, VariantCall
from resolwe_bio.variants.views import VariantAnnotationSerializer


class AnnotationTestCase(BioProcessTestCase):
    @tag_process("variant-annotation")
    def test_variant_annotation(self):
        with self.preparation_stage():
            self.contributor = User.objects.get_or_create(
                username="contributor_variants",
                email="contributor_variants@genialis.com",
            )[0]

            self.sample = Sample.objects.create(contributor=self.contributor)

            input_folder = Path("variant_annotation") / "input"
            output_folder = Path("variant_annotation") / "output"
            dbsnp = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "dbsnp.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            clinvar = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "clinvar.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            canonical_transcripts = self.run_process(
                "upload-geneset",
                {
                    "src": input_folder / "canonical_transcripts_ens109.txt.gz",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                },
            )
            dbsnp_2 = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "dbsnp_2.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )
            clinvar_2 = self.run_process(
                "upload-variants-vcf",
                {
                    "src": input_folder / "clinvar_2.vcf.gz",
                    "species": "Homo sapiens",
                    "build": "GRCh38",
                },
            )

        variant_annotation = self.run_process(
            "variant-annotation",
            {
                "database": "GRCh38.109",
                "dbsnp": dbsnp.id,
                "clinvar": clinvar.id,
                "canonical_transcripts": canonical_transcripts.id,
                "genome_assembly": "GRCh38",
            },
            Data.STATUS_ERROR,
        )

        error_msg = ["No variants found in the database for the selected inputs."]
        self.assertEqual(variant_annotation.process_error, error_msg)

        variant_T_A = Variant.objects.get_or_create(
            species="Homo sapiens",
            genome_assembly="GRCh38",
            chromosome="13",
            position=48307307,
            reference="T",
            alternative="A",
        )
        VariantCall.objects.get_or_create(
            sample=self.sample,
            variant=variant_T_A[0],
            quality=70,
            depth_norm_quality=0.7,
            alternative_allele_depth=11,
            depth=30,
            filter="PASS",
            genotype="0/1",
            genotype_quality=99,
        )

        variant_A_T = Variant.objects.get_or_create(
            species="Homo sapiens",
            genome_assembly="GRCh38",
            chromosome="12",
            position=25206394,
            reference="A",
            alternative="T",
        )
        VariantCall.objects.get_or_create(
            sample=self.sample,
            variant=variant_A_T[0],
            quality=70,
            depth_norm_quality=0.7,
            alternative_allele_depth=11,
            depth=30,
            filter="PASS",
            genotype="0/1",
            genotype_quality=99,
        )

        variant_annotation = self.run_process(
            "variant-annotation",
            {
                "database": "GRCh38.109",
                "dbsnp": dbsnp.id,
                "clinvar": clinvar.id,
                "canonical_transcripts": canonical_transcripts.id,
                "genome_assembly": "GRCh38",
                "advanced": {"save_vcf": True},
            },
        )

        expected_transcript = {
            "id": variant_A_T[0].annotation.transcripts.first().id,
            "variant_annotation_id": variant_A_T[0].annotation.pk,
            "annotation": "3_prime_UTR_variant",
            "annotation_impact": "MODIFIER",
            "gene": "KRAS",
            "protein_impact": "",
            "transcript_id": "ENST00000256078",
            "canonical": False,
        }

        expected_annotation = {
            "id": variant_A_T[0].annotation.id,
            "variant_id": variant_A_T[0].id,
            "type": "SNP",
            "clinical_diagnosis": "Noonan_syndrome",
            "clinical_significance": "Likely_benign",
            "dbsnp_id": "rs1137189",
            "clinvar_id": "308075",
        }

        serialized_annotation = VariantAnnotationSerializer(
            variant_A_T[0].annotation
        ).data
        self.assertDictEqual(
            serialized_annotation["transcripts"][0],
            expected_transcript,
        )

        del serialized_annotation["transcripts"]
        self.assertDictEqual(
            serialized_annotation,
            expected_annotation,
        )

        self.assertFile(
            variant_annotation, "vcf", output_folder / "annotated_variants.vcf.gz"
        )

        # This test is added to cover the case when the user want to override the existing annotation.
        variant_annotation = self.run_process(
            "variant-annotation",
            {
                "database": "GRCh38.109",
                "dbsnp": dbsnp_2.id,
                "clinvar": clinvar_2.id,
                "canonical_transcripts": canonical_transcripts.id,
                "genome_assembly": "GRCh38",
                "select_all": True,
                "advanced": {"save_vcf": True},
            },
        )

        variant_A_T[0].refresh_from_db()

        serialized_annotation = VariantAnnotationSerializer(
            variant_A_T[0].annotation
        ).data

        expected_annotation["id"] = serialized_annotation["id"]
        expected_annotation["dbsnp_id"] = ""
        expected_annotation["clinvar_id"] = ""

        del serialized_annotation["transcripts"]

        self.assertDictEqual(
            serialized_annotation,
            expected_annotation,
        )
