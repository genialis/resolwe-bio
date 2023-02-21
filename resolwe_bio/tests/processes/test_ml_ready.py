from pathlib import Path

from resolwe.flow.models import Collection
from resolwe.test import tag_process, with_resolwe_host

from resolwe_bio.models import Sample
from resolwe_bio.utils.test import TEST_FILES_DIR, KBBioProcessTestCase


class MLReadyTestCase(KBBioProcessTestCase):
    def prepare_collection_with_one_sample(self):
        """Prepare a collection with one sample."""
        collection = Collection.objects.create(
            name="Reference spaces of Data Flywheel for AI/ML",
            slug="reference-spaces",
            contributor=self.contributor,
        )

        sample_1 = Sample.objects.create(name="Sample 1", contributor=self.contributor)
        sample_1.slug = "sample-1"
        sample_1.name = "SAMPLE 1"
        sample_1.collection = collection
        sample_1.save()

        return collection, sample_1

    def prepare_exp_matrix(self, sample_ids, feature_ids, outfile):
        """
        Prepare exp matrix.

        This matrix needs to be prepared on the fly, as it contains
        sample IDs that are determined on the fly as well.
        """
        # Create matrix
        exp_matrix = [["sample_id"] + list(feature_ids)]
        for sample_id in sample_ids:
            exp_matrix.append(
                [sample_id] + [sample_id * i for i, _ in enumerate(feature_ids)]
            )

        # Write to file
        with open(outfile, "wt") as handle:
            for row in exp_matrix:
                handle.write("\t".join(map(str, row)) + "\n")

    @with_resolwe_host
    @tag_process("reference-space")
    def test_upload_reference_space(self):
        base = Path("ml_ready")
        inputs = base / "input"
        outputs = base / "output"
        with self.preparation_stage():
            collection, sample_1 = self.prepare_collection_with_one_sample()

            self.prepare_exp_matrix(
                sample_ids=[sample_1.id],
                feature_ids=["ENSG00000178591", "ENSG00000101255"],
                outfile=str(TEST_FILES_DIR / inputs / "exp_matrix.tsv"),
            )

        ref_space = self.run_process(
            "reference-space",
            {
                "name": "My name",
                "description": "My description",
                "source": "ENSEMBL",
                "species": "Homo sapiens",
                "training_data": str(inputs / "exp_matrix.tsv"),
                "preprocessor": str(inputs / "preprocessor.pickle"),
            },
            collection=collection,
        )

        self.assertEqual(ref_space.output["source"], "ENSEMBL")
        self.assertEqual(ref_space.output["species"], "Homo sapiens")
        self.assertFile(ref_space, "preprocessor", str(inputs / "preprocessor.pickle"))
        self.assertJSON(
            ref_space,
            ref_space.output["features"],
            "",
            str(outputs / "features.json.gz"),
        )
        self.assertEqual(ref_space.descriptor, {"description": "My description"})
        self.assertEqual(ref_space.descriptor_schema.slug, "geneset")

    @tag_process("reference-space", "upload-ml-expression")
    def test_upload_ml_expression(self):
        base = Path("ml_ready")
        inputs = base / "input"

        with self.preparation_stage():
            collection, sample_1 = self.prepare_collection_with_one_sample()

            self.prepare_exp_matrix(
                sample_ids=[sample_1.id],
                feature_ids=["ENSG00000178591", "ENSG00000101255"],
                outfile=str(TEST_FILES_DIR / inputs / "exp_matrix.tsv"),
            )

            ref_space = self.run_process(
                "reference-space",
                {
                    "name": "My name",
                    "description": "My description",
                    "source": "ENSEMBL",
                    "species": "Homo sapiens",
                    "training_data": str(inputs / "exp_matrix.tsv"),
                    "preprocessor": str(inputs / "preprocessor.pickle"),
                },
                collection=collection,
            )

        ml_exp = self.run_process(
            "upload-ml-expression",
            {
                "source": "ENSEMBL",
                "species": "Homo sapiens",
                "exp": str(inputs / "exp_matrix.tsv"),
                "reference_space": ref_space.id,
            },
            collection=collection,
        )

        self.assertEqual(ml_exp.output["source"], "ENSEMBL")
        self.assertEqual(ml_exp.output["species"], "Homo sapiens")
        # Cannot test with assertFile, like:
        # self.assertFile(ml_exp, "exp", outputs / "exp_matrix.tsv")
        # Since it is not possible to predict sample ID during tests

        # Get the original data
        with open(str(TEST_FILES_DIR / inputs / "exp_matrix.tsv")) as handle:
            contents_original = [line.strip().split("\t") for line in handle]

        # Get the processed data
        exp_path = ml_exp.location.get_path(filename=ml_exp.output["exp"]["file"])
        with open(exp_path) as handle:
            contents_test = [line.strip().split("\t") for line in handle]

        self.assertTrue(len(contents_test) == len(contents_original) == 2)
        # Header needs to start with sample_id and than sorted features
        self.assertEqual(
            contents_test[0],
            [contents_original[0][0]] + sorted(contents_original[0][1:]),
        )
        # Only one line, test except for sample id
        self.assertCountEqual(contents_test[1][1:], contents_original[1][1:])
