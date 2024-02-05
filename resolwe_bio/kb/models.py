""".. Ignore pydocstyle D400.

======
Models
======

"""

from django.contrib.postgres.fields import ArrayField
from django.contrib.postgres.indexes import GinIndex
from django.contrib.postgres.search import SearchVectorField
from django.db import models


# NOTE: Feature are manually inserted in management command, so take
#       care that it is synced with model definition.
class Feature(models.Model):
    """Describes a feature in the knowledge base."""

    TYPE_GENE = "gene"
    TYPE_TRANSCRIPT = "transcript"
    TYPE_EXON = "exon"
    TYPE_PROBE = "probe"
    TYPE_CHOICES = (
        (TYPE_GENE, "Gene"),
        (TYPE_TRANSCRIPT, "Transcript"),
        (TYPE_EXON, "Exon"),
        (TYPE_PROBE, "Probe"),
    )

    SUBTYPE_PROTEIN_CODING = "protein-coding"
    SUBTYPE_PSEUDO = "pseudo"
    SUBTYPE_RRNA = "rRNA"
    SUBTYPE_NCRNA = "ncRNA"
    SUBTYPE_SNRNA = "snRNA"
    SUBTYPE_SNORNA = "snoRNA"
    SUBTYPE_TRNA = "tRNA"
    SUBTYPE_ASRNA = "asRNA"
    SUBTYPE_OTHER = "other"
    SUBTYPE_UNKNOWN = "unknown"
    SUBTYPE_CHOICES = (
        (SUBTYPE_PROTEIN_CODING, "Protein-coding"),
        (SUBTYPE_PSEUDO, "Pseudo"),
        (SUBTYPE_RRNA, "rRNA"),
        (SUBTYPE_NCRNA, "ncRNA"),
        (SUBTYPE_SNRNA, "snRNA"),
        (SUBTYPE_SNORNA, "snoRNA"),
        (SUBTYPE_TRNA, "tRNA"),
        (SUBTYPE_ASRNA, "asRNA"),
        (SUBTYPE_OTHER, "Other"),
        (SUBTYPE_UNKNOWN, "Unknown"),
    )

    # Because Django ORM cannot handle composite primary keys, each feature is
    # still assigned an internal numeric 'id' and the ('source', 'feature_id',
    # 'species') combination is used to uniquely identify a feature.
    source = models.CharField(max_length=20)
    feature_id = models.CharField(max_length=50)
    species = models.CharField(max_length=50)

    type = models.CharField(max_length=20, choices=TYPE_CHOICES)
    sub_type = models.CharField(max_length=20, choices=SUBTYPE_CHOICES)
    name = models.CharField(max_length=1024)
    full_name = models.CharField(max_length=350, blank=True)
    description = models.TextField(blank=True)
    aliases = ArrayField(models.CharField(max_length=256), default=list, blank=True)

    #: field used for full-text search
    search = SearchVectorField(null=True)

    class Meta:
        """Feature Meta options."""

        constraints = [
            models.UniqueConstraint(
                fields=["source", "feature_id", "species"],
                name="uniq_feature_source_feature_id_species",
            ),
        ]
        indexes = [
            models.Index(name="idx_feature_source", fields=["source"]),
            models.Index(name="idx_feature_species", fields=["species"]),
            models.Index(name="idx_feature_feature_id", fields=["feature_id"]),
            models.Index(name="idx_feature_type", fields=["type"]),
            GinIndex(name="idx_feature_search", fields=["search"]),
        ]

    def __str__(self):
        """Represent a feature instance as a string."""
        return "{source}: {feature_id} ({species})".format(
            source=self.source,
            feature_id=self.feature_id,
            species=self.species,
        )


# NOTE: Mappings are manually inserted in management command, so take
#       care that it is synced with model definition.
class Mapping(models.Model):
    """Describes a mapping between features from different sources."""

    RELATION_TYPE_CROSSDB = "crossdb"
    RELATION_TYPE_ORTHOLOG = "ortholog"
    RELATION_TYPE_TRANSCRIPT = "transcript"
    RELATION_TYPE_EXON = "exon"
    RELATION_TYPE_CHOICES = (
        (RELATION_TYPE_CROSSDB, "Crossdb"),
        (RELATION_TYPE_ORTHOLOG, "Ortholog"),
        (RELATION_TYPE_TRANSCRIPT, "Transcript"),
        (RELATION_TYPE_EXON, "Exon"),
    )

    relation_type = models.CharField(max_length=20, choices=RELATION_TYPE_CHOICES)
    source_db = models.CharField(max_length=20)
    source_id = models.CharField(max_length=50)
    source_species = models.CharField(max_length=50)
    target_db = models.CharField(max_length=20)
    target_id = models.CharField(max_length=50)
    target_species = models.CharField(max_length=50)

    class Meta:
        """Mapping Meta options."""

        constraints = [
            models.UniqueConstraint(
                fields=[
                    "source_db",
                    "source_id",
                    "source_species",
                    "target_db",
                    "target_id",
                    "target_species",
                    "relation_type",
                ],
                name="uniq_mapping_source_target_type",
            ),
        ]
        indexes = [
            models.Index(
                name="idx_feature_source_target",
                fields=[
                    "source_db",
                    "source_id",
                    "source_species",
                    "target_db",
                    "target_species",
                ],
            ),
            models.Index(
                name="idx_feature_target",
                fields=["target_db", "target_id", "target_species"],
            ),
        ]

    def __str__(self):
        """Represent a feature instance as a string."""
        return "{src_db}: {src_id} ({src_species}) -> {dst_db}: {dst_id} ({dst_species})".format(
            src_db=self.source_db,
            src_id=self.source_id,
            src_species=self.source_species,
            dst_db=self.target_db,
            dst_id=self.target_id,
            dst_species=self.target_species,
        )
