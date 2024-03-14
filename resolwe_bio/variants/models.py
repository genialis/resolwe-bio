""".. Ignore pydocstyle D400.

======
Models
======

"""

from django.conf import settings
from django.contrib.postgres.fields import ArrayField
from django.db import models

from resolwe.flow.models import Data
from resolwe.flow.models import Entity as Sample

# TODO: sync with the kb, at least the first entry.
SPECIES_MAX_LENGTH = 50
GENOME_ASSEMBLY_MAX_LENGTH = 20
CHROMOSOME_MAX_LENGTH = 20
REFERENCE_MAX_LENGTH = 100
ALTERNATIVE_MAX_LENGTH = 100

# Metadata
VARIANT_DATA_SOURCE_MAX_LENGTH = 100
VARIANT_ANNOTATION_SOURCE_MAX_LENGTH = 100
ANNOTATION_MAX_LENGTH = 200
GENOTYPE_MAX_LENGTH = 100
TYPE_MAX_LENGTH = 100
CLINICAL_DIAGNOSIS_MAX_LENGTH = 200
CLINICAL_SIGNIFICANCE_MAX_LENGTH = 100
DBSNP_ID_MAX_LENGTH = 20
CLINICAL_VAR_ID_MAX_LENGTH = 20
ANNOTATION_IMPACT_MAX_LENGTH = 20
GENE_MAX_LENGTH = 100
PROTEIN_IMPACT_MAX_LENGTH = 100
FEATURE_ID_MAX_LENGTH = 200
FILTER_MAX_LENGTH = 20


class Variant(models.Model):
    """Describe a variant in the database."""

    # Because Django ORM cannot handle composite primary keys, each feature is
    # still assigned an internal numeric 'id' and the ('source', 'feature_id',
    # 'species') combination is used to uniquely identify a feature.
    species = models.CharField(max_length=SPECIES_MAX_LENGTH)
    genome_assembly = models.CharField(max_length=GENOME_ASSEMBLY_MAX_LENGTH)
    chromosome = models.CharField(max_length=CHROMOSOME_MAX_LENGTH)
    position = models.PositiveBigIntegerField()
    reference = models.CharField(max_length=REFERENCE_MAX_LENGTH)
    alternative = models.CharField(max_length=ALTERNATIVE_MAX_LENGTH)

    class Meta:
        """Add constraint for composite key."""

        constraints = [
            models.UniqueConstraint(
                fields=[
                    "species",
                    "genome_assembly",
                    "chromosome",
                    "position",
                    "reference",
                    "alternative",
                ],
                name="uniq_composite_key_variants",
            ),
        ]


class VariantAnnotation(models.Model):
    """Describes an annotation of a variant."""

    variant = models.OneToOneField(
        Variant, on_delete=models.CASCADE, related_name="annotation"
    )
    type = models.CharField(max_length=TYPE_MAX_LENGTH, blank=True, null=True)
    annotation = models.CharField(max_length=ANNOTATION_MAX_LENGTH)
    annotation_impact = models.CharField(max_length=ANNOTATION_IMPACT_MAX_LENGTH)
    gene = models.CharField(max_length=GENE_MAX_LENGTH)
    protein_impact = models.CharField(max_length=PROTEIN_IMPACT_MAX_LENGTH)
    feature_id = ArrayField(
        models.CharField(max_length=FEATURE_ID_MAX_LENGTH), default=list
    )
    clinical_diagnosis = models.CharField(
        max_length=CLINICAL_DIAGNOSIS_MAX_LENGTH, blank=True, null=True
    )
    clinical_significance = models.CharField(
        max_length=CLINICAL_SIGNIFICANCE_MAX_LENGTH, blank=True, null=True
    )
    # Optional references to the external databases.
    dbsnp_id = models.CharField(max_length=DBSNP_ID_MAX_LENGTH, blank=True, null=True)
    clinical_var_id = models.CharField(
        max_length=CLINICAL_VAR_ID_MAX_LENGTH, blank=True, null=True
    )
    data = models.ForeignKey(
        Data,
        on_delete=models.CASCADE,
        related_name="variant_annotations",
        null=True,
        blank=True,
    )


class VariantExperiment(models.Model):
    """Represents a single experiment."""

    variant_data_source = models.CharField(max_length=VARIANT_DATA_SOURCE_MAX_LENGTH)
    date = models.DateTimeField(auto_now_add=True, db_index=True)
    contributor = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.PROTECT)


class VariantCall(models.Model):
    """VariantCall object."""

    sample = models.ForeignKey(
        Sample,
        on_delete=models.CASCADE,
        related_name="variant_calls",
        null=True,
        blank=True,
    )

    variant = models.ForeignKey(
        Variant,
        on_delete=models.CASCADE,
        related_name="variant_calls",
        null=True,
        blank=True,
    )

    experiment = models.ForeignKey(
        VariantExperiment,
        on_delete=models.CASCADE,
        related_name="variant_calls",
        null=True,
        blank=True,
    )

    # QC data.
    quality = models.FloatField()
    depth = models.PositiveIntegerField()
    filter = models.CharField(max_length=FILTER_MAX_LENGTH)
    genotype = models.CharField(max_length=GENOTYPE_MAX_LENGTH, blank=True, null=True)
    data = models.ForeignKey(
        Data,
        on_delete=models.CASCADE,
        related_name="variant_calls",
        null=True,
        blank=True,
    )
