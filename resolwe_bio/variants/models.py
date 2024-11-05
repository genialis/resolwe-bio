""".. Ignore pydocstyle D400.

======
Models
======

"""

from django.conf import settings
from django.db import models

from resolwe.auditlog.models import AuditModel
from resolwe.flow.models import Data
from resolwe.flow.models import Entity as Sample
from resolwe.permissions.models import PermissionInterface

# TODO: sync with the kb, at least the first entry.
SPECIES_MAX_LENGTH = 50
GENOME_ASSEMBLY_MAX_LENGTH = 20
CHROMOSOME_MAX_LENGTH = 20

# Metadata
VARIANT_DATA_SOURCE_MAX_LENGTH = 100
VARIANT_ANNOTATION_SOURCE_MAX_LENGTH = 100
GENOTYPE_MAX_LENGTH = 3
TYPE_MAX_LENGTH = 100
CLINICAL_SIGNIFICANCE_MAX_LENGTH = 100
DBSNP_ID_MAX_LENGTH = 150
CLINICAL_VAR_ID_MAX_LENGTH = 150
ANNOTATION_IMPACT_MAX_LENGTH = 20
GENE_MAX_LENGTH = 100
PROTEIN_IMPACT_MAX_LENGTH = 100
FEATURE_ID_MAX_LENGTH = 200
FILTER_MAX_LENGTH = 20

VariantPrimaryKey = (
    "species",
    "genome_assembly",
    "chromosome",
    "position",
    "reference",
    "alternative",
)


class Variant(AuditModel):
    """Describe a variant in the database."""

    # Because Django ORM cannot handle composite primary keys, each feature is
    # still assigned an internal numeric 'id' and the ('species', 'genome_assembly',
    # 'chromosome', 'position', 'reference', 'alternative) combination is used to
    # uniquely identify a variant.

    #: species name
    species = models.CharField(max_length=SPECIES_MAX_LENGTH)

    #: genome assembly
    genome_assembly = models.CharField(max_length=GENOME_ASSEMBLY_MAX_LENGTH)

    #: chromosome
    chromosome = models.CharField(max_length=CHROMOSOME_MAX_LENGTH)

    #: position
    position = models.PositiveBigIntegerField()

    #: reference
    reference = models.TextField()

    #: alternative
    alternative = models.TextField()

    class Meta:
        """Add constraint for composite key."""

        constraints = [
            models.UniqueConstraint(
                fields=VariantPrimaryKey,
                name="uniq_composite_key_variants",
            ),
        ]


class VariantAnnotation(AuditModel):
    """Describes an annotation of a variant."""

    TYPE_CHOICES = (
        ("SNP", "SNP"),
        ("INDEL", "INDEL"),
        ("MIXED", "MIXED"),
        ("MNP", "MNP"),
        ("SYMBOLIC", "SYMBOLIC"),
    )

    #: the referenced variant
    variant = models.OneToOneField(
        Variant, on_delete=models.CASCADE, related_name="annotation"
    )

    #: annotation type
    type = models.CharField(
        max_length=TYPE_MAX_LENGTH, blank=True, null=True, choices=TYPE_CHOICES
    )

    #: clinical diagnosis
    clinical_diagnosis = models.TextField(blank=True, null=True)

    #: clinical significance
    clinical_significance = models.CharField(
        max_length=CLINICAL_SIGNIFICANCE_MAX_LENGTH, blank=True, null=True
    )

    #: reference to dbsnp external database
    dbsnp_id = models.CharField(max_length=DBSNP_ID_MAX_LENGTH, blank=True, null=True)

    #: reference to clinvar external database
    clinvar_id = models.CharField(
        max_length=CLINICAL_VAR_ID_MAX_LENGTH, blank=True, null=True
    )

    #: reference to the data object
    data = models.ForeignKey(
        Data,
        on_delete=models.CASCADE,
        related_name="variant_annotations",
        null=True,
        blank=True,
    )


class VariantAnnotationTranscript(AuditModel):
    """Represents a transcript in a variant annotation."""

    #: referenced variant annotation model
    variant_annotation = models.ForeignKey(
        "VariantAnnotation", on_delete=models.CASCADE, related_name="transcripts"
    )

    #: annotation
    annotation = models.TextField()

    #: impact
    annotation_impact = models.CharField(max_length=ANNOTATION_IMPACT_MAX_LENGTH)

    #: gene
    gene = models.CharField(max_length=GENE_MAX_LENGTH)

    #: protein impact
    protein_impact = models.CharField(max_length=PROTEIN_IMPACT_MAX_LENGTH)

    #: transcript_id
    transcript_id = models.CharField(max_length=FEATURE_ID_MAX_LENGTH)

    #: is this transcript canonical
    canonical = models.BooleanField(default=False)


class VariantExperiment(AuditModel):
    """Represents a single experiment."""

    #: data source
    variant_data_source = models.CharField(max_length=VARIANT_DATA_SOURCE_MAX_LENGTH)

    #: timestamp of the experiment
    timestamp = models.DateTimeField(auto_now_add=True, db_index=True)

    #: source of the annotation
    contributor = models.ForeignKey(settings.AUTH_USER_MODEL, on_delete=models.PROTECT)


class VariantCall(AuditModel, PermissionInterface):
    """VariantCall object."""

    @classmethod
    def permission_proxy(cls) -> str:
        """Return the permission proxy name."""
        return "sample"

    #: the referenced sample
    sample = models.ForeignKey(
        Sample, on_delete=models.CASCADE, related_name="variant_calls"
    )

    #: the referenced variant
    variant = models.ForeignKey(
        Variant,
        on_delete=models.CASCADE,
        related_name="variant_calls",
        null=True,
        blank=True,
    )

    #: the referenced experiment
    experiment = models.ForeignKey(
        VariantExperiment,
        on_delete=models.CASCADE,
        related_name="variant_calls",
        null=True,
        blank=True,
    )

    #: quality
    quality = models.FloatField(blank=True, null=True)

    #: depth_norm_quality
    depth_norm_quality = models.FloatField(blank=True, null=True)

    #: alternative_allele_depth
    alternative_allele_depth = models.PositiveIntegerField(blank=True, null=True)

    #: depth
    depth = models.PositiveIntegerField(blank=True, null=True)

    #: genotype
    genotype = models.CharField(max_length=GENOTYPE_MAX_LENGTH, blank=True, null=True)

    #: genotype quality
    genotype_quality = models.PositiveIntegerField(blank=True, null=True)

    #: filter
    filter = models.CharField(max_length=FILTER_MAX_LENGTH, blank=True, null=True)

    #: the referenced data object
    data = models.ForeignKey(
        Data,
        on_delete=models.CASCADE,
        related_name="variant_calls",
        null=True,
        blank=True,
    )
