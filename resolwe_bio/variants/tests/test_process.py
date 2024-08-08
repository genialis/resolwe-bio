from resolwe.process import IntegerField, ListField, StringField

from resolwe_bio.process.runtime import ProcessBio


class SimpleGetVariant(ProcessBio):
    slug = "test-variants-simple-get"
    name = "Test variants"
    process_type = "data:varianttest"
    version = "1.0.0"

    class Output:
        """Output fields."""

        variant_species = ListField(StringField(), label="Variants", required=True)

    def run(self, inputs, outputs):
        """Start the process."""
        variant_species = [variant.species for variant in self.variant.iterate()]
        outputs.variant_species = variant_species


class SimpleFilterVariant(ProcessBio):
    slug = "test-variants-simple-filter"
    name = "Test variants"
    process_type = "data:varianttest"
    version = "1.0.0"

    class Output:
        """Output fields."""

        variant_species = ListField(StringField(), label="Variants", required=True)

    class Input:
        """Input fields."""

        species_filter = StringField(label="Species filter")

    def run(self, inputs, outputs):
        """Start the process."""
        variant_species = [
            variant.species
            for variant in self.variant.iterate(species=inputs.species_filter)
        ]
        outputs.variant_species = variant_species


class SimpleGetVariantCall(ProcessBio):
    slug = "test-variantcalls-simple-get"
    name = "Test variants"
    process_type = "data:varianttest"
    version = "1.0.0"

    class Output:
        """Output fields."""

        filters = ListField(StringField(), label="Filters", required=True)

    def run(self, inputs, outputs):
        """Start the process."""
        filters = [variant_call.filter for variant_call in self.variant_call.iterate()]
        outputs.filters = filters


class GetVariantPosition(ProcessBio):
    slug = "test-variant-position-get"
    name = "Test variants"
    process_type = "data:varianttest"
    version = "1.0.0"

    class Output:
        """Output fields."""

        position = IntegerField(label="Position", required=True)

    class Input:
        """Input fields."""

        variant_id = IntegerField(label="Variant id")

    def run(self, inputs, outputs):
        """Start the process."""
        outputs.position = self.variant.get(id=inputs.variant_id).position
