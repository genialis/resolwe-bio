"""Upload methylation array data (IDAT)."""
from resolwe.process import FileField, Process, SchedulingClass, StringField


def validate_filename_suffix(filename, suffix, resolwe_process=Process):
    """Raise an error if unexpected file name suffix is encountered."""
    try:
        assert filename.endswith(suffix)
    except AssertionError:
        resolwe_process.error(
            f"Unsupported file name extension. A file {filename} "
            f"should end with {suffix}."
        )


class UploadIdatData(Process):
    """Upload Illumina methylation array raw IDAT data.

    This import process accepts Illumina methylation array BeadChip raw
    files in IDAT format. Two input files, one for each of the Green and
    Red signal channels, are expected. The uploads of human (HM27, HM450,
    EPIC) and mouse (MM285) array types are supported.
    """

    slug = "upload-idat"
    name = "IDAT file"
    process_type = "data:methylationarray:idat"
    version = "1.0.0"
    category = "Import"
    scheduling_class = SchedulingClass.BATCH
    requirements = {
        "expression-engine": "jinja",
        "executor": {
            "docker": {"image": "public.ecr.aws/s4q6j6e8/resolwebio/common:2.7.0"}
        },
        "resources": {"cores": 1, "memory": 2048},
    }
    entity = {
        "type": "sample",
        "descriptor_schema": "sample",
    }
    data_name = "{{ red_channel.file|default('?') }}"

    class Input:
        """Input field to process UploadIdatData."""

        red_channel = FileField(label="Red channel IDAT file (*_Red.idat)")
        green_channel = FileField(label="Green channel IDAT file (*_Grn.idat)")
        species = StringField(
            label="Species",
            description="Select a species name from the dropdown menu.",
            default="Homo sapiens",
            choices=[
                ("Homo sapiens", "Homo sapiens"),
                ("Mus musculus", "Mus musculus"),
            ],
        )
        platform = StringField(
            label="Protein ID database source",
            description="Select a methylation array platform for human "
            "(HM450, HM27, EPIC) or mouse (MM285) samples.",
            default="HM450",
            choices=[
                ("HM450", "HM450"),
                ("HM27", "HM27"),
                ("EPIC", "EPIC"),
                ("MM285", "MM285"),
            ],
        )

    class Output:
        """Output field of the process UploadProteomicsData."""

        red_channel = FileField(label="Red channel IDAT file")
        green_channel = FileField(label="Green channel IDAT file")
        species = StringField(label="Species")
        platform = StringField(label="Platform")

    def run(self, inputs, outputs):
        """Run the analysis."""

        if inputs.species == "Mus musculus" and inputs.platform != "MM285":
            self.error(
                f"Platform type {inputs.platform} does not match the selected species {inputs.species}."
            )

        red = inputs.red_channel.import_file(imported_format="compressed")
        grn = inputs.green_channel.import_file(imported_format="compressed")

        validate_filename_suffix(red, "_Red.idat.gz")
        validate_filename_suffix(grn, "_Grn.idat.gz")

        sample_name_red = red[:-12]
        sample_name_grn = grn[:-12]

        if sample_name_red != sample_name_grn:
            self.error(
                "The input IDAT files don't have a matching filename prefix. "
                "The sample data might be mismatched."
            )

        outputs.red_channel = red
        outputs.green_channel = grn
        outputs.species = inputs.species
        outputs.platform = inputs.platform
