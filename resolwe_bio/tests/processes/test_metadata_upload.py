from pathlib import Path

from resolwe.flow.models import Collection, Data
from resolwe.test import tag_process

from resolwe_bio.models import Sample
from resolwe_bio.utils.test import BioProcessTestCase


class MetaDataUploadProcessorTestCase(BioProcessTestCase):
    @tag_process("upload-metadata-unique")
    def test_upload_metadata_unique(self):
        base = Path("metadata_upload")
        inputs = base / "inputs"

        collection = Collection.objects.create(
            name="Test collection", contributor=self.contributor
        )

        sample_1 = Sample.objects.create(name="Sample 1", contributor=self.contributor)
        sample_1.slug = "sample-1"
        sample_1.name = "SAMPLE 1"
        sample_1.collection = collection
        sample_1.save()

        sample_2 = Sample.objects.create(name="Sample 2", contributor=self.contributor)
        sample_2.slug = "sample-2"
        sample_2.name = "SAMPLE 2"
        sample_2.collection = collection
        sample_2.save()

        meta = self.run_process(
            "upload-metadata-unique",
            {
                "src": str(inputs / "sample_metatable.tsv"),
            },
            collection=collection,
        )

        self.assertFile(meta, "table", str(inputs / "sample_metatable.tsv"))
        self.assertFields(meta, "n_samples", 8)

        upper_case_suffix = self.run_process(
            "upload-metadata-unique",
            {
                "src": str(inputs / "upper_case_sufix.TSV"),
            },
            collection=collection,
        )

        info = ["File extension of the table was replaced with a lower case version."]
        self.assertEqual(upper_case_suffix.process_info, info)

        mixed_case_suffix = self.run_process(
            "upload-metadata-unique",
            {
                "src": str(inputs / "mixed_case_sufix.TsV"),
            },
            collection=collection,
        )

        info = ["File extension of the table was replaced with a lower case version."]
        self.assertEqual(mixed_case_suffix.process_info, info)

        meta_xlsx = self.run_process(
            "upload-metadata-unique",
            {
                "src": str(inputs / "sample_metatable.xlsx"),
            },
            collection=collection,
        )
        self.assertFile(meta_xlsx, "table", str(inputs / "sample_metatable.xlsx"))
        self.assertFields(meta_xlsx, "n_samples", 8)

        three_header_format = self.run_process(
            "upload-metadata-unique",
            {
                "src": str(inputs / "iris_legacy.TAB"),
            },
            Data.STATUS_ERROR,
            collection=collection,
        )
        error_msg = [
            "The uploaded metadata table needs to contain "
            "exactly one of the following columns: "
            "['Sample ID', 'Sample name', 'Sample slug']."
        ]
        self.assertEqual(three_header_format.process_error, error_msg)

        wo_collection = self.run_process(
            "upload-metadata-unique",
            {
                "src": str(inputs / "sample_metatable.tsv"),
            },
            Data.STATUS_ERROR,
        )
        error_msg = [
            "Metadata table was not uploaded to a Collection. "
            "Matching of metadata entries to Sample objects is not possible."
        ]
        self.assertEqual(wo_collection.process_error, error_msg)

        one_to_many = self.run_process(
            "upload-metadata-unique",
            {
                "src": str(inputs / "sample_metatable_one_to_many.tsv"),
            },
            Data.STATUS_ERROR,
            collection=collection,
        )
        error_msg = [
            "Duplicated metadata entries ['sample-1'] were found. "
            "Please use the metadata upload process that "
            "allows for one-to-many relations instead."
        ]
        self.assertEqual(one_to_many.process_error, error_msg)

        empty = self.run_process(
            "upload-metadata-unique",
            {
                "src": str(inputs / "empty.tsv"),
            },
            Data.STATUS_ERROR,
            collection=collection,
        )
        error_msg = [
            "It was not possible to read the provided data table. "
            "No columns to parse from file"
        ]
        self.assertEqual(empty.process_error, error_msg)

    @tag_process("upload-metadata")
    def test_upload_metadata(self):
        base = Path("metadata_upload")
        inputs = base / "inputs"

        collection = Collection.objects.create(
            name="Test collection", contributor=self.contributor
        )

        sample_1 = Sample.objects.create(name="Sample 1", contributor=self.contributor)
        sample_1.slug = "sample-1"
        sample_1.name = "SAMPLE 1"
        sample_1.collection = collection
        sample_1.save()

        sample_2 = Sample.objects.create(name="Sample 2", contributor=self.contributor)
        sample_2.slug = "sample-2"
        sample_2.name = "SAMPLE 2"
        sample_2.collection = collection
        sample_2.save()

        one_to_many = self.run_process(
            "upload-metadata",
            {
                "src": str(inputs / "sample_metatable_one_to_many.tsv"),
            },
            collection=collection,
        )

        self.assertFile(
            one_to_many, "table", str(inputs / "sample_metatable_one_to_many.tsv")
        )
        self.assertFields(one_to_many, "n_samples", 7)
