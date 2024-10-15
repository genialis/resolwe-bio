""".. Ignore pydocstyle D400.

===============
Variant Filters
===============

"""

import django_filters as filters

from resolwe.flow.filters import (
    DATETIME_LOOKUPS,
    NUMBER_LOOKUPS,
    RELATED_LOOKUPS,
    TEXT_LOOKUPS,
    CheckQueryParamsMixin,
)

from resolwe_bio.variants.models import (
    Variant,
    VariantAnnotation,
    VariantCall,
    VariantExperiment,
)


class VariantFilter(CheckQueryParamsMixin, filters.FilterSet):
    """Filter the Variant objects endpoint."""

    class Meta:
        """Filter configuration."""

        model = Variant
        fields = {
            "id": NUMBER_LOOKUPS,
            "species": TEXT_LOOKUPS,
            "genome_assembly": TEXT_LOOKUPS,
            "chromosome": TEXT_LOOKUPS,
            "position": NUMBER_LOOKUPS,
            "reference": TEXT_LOOKUPS,
            "alternative": TEXT_LOOKUPS,
            "annotation__type": TEXT_LOOKUPS,
            "annotation__transcripts": RELATED_LOOKUPS,
            "annotation__transcripts__annotation": TEXT_LOOKUPS,
            "annotation__transcripts__annotation_impact": TEXT_LOOKUPS,
            "annotation__transcripts__gene": TEXT_LOOKUPS,
            "annotation__transcripts__protein_impact": TEXT_LOOKUPS,
            "annotation__transcripts__canonical": ["exact"],
            "annotation__type": TEXT_LOOKUPS,
            "annotation__clinical_diagnosis": TEXT_LOOKUPS,
            "annotation__clinical_significance": TEXT_LOOKUPS,
            "annotation__dbsnp_id": TEXT_LOOKUPS,
            "annotation__clinvar_id": TEXT_LOOKUPS,
            "annotation__data": RELATED_LOOKUPS,
            "variant_calls": RELATED_LOOKUPS,
            "variant_calls__quality": NUMBER_LOOKUPS,
            "variant_calls__depth": NUMBER_LOOKUPS,
            "variant_calls__sample__slug": TEXT_LOOKUPS,
            "variant_calls__sample": RELATED_LOOKUPS,
            "variant_calls__experiment": RELATED_LOOKUPS,
            "variant_calls__filter": TEXT_LOOKUPS,
            "variant_calls__genotype": TEXT_LOOKUPS,
            "variant_calls__genotype_quality": NUMBER_LOOKUPS,
            "variant_calls__alternative_allele_depth": NUMBER_LOOKUPS,
            "variant_calls__depth_norm_quality": NUMBER_LOOKUPS,
        }

    @property
    def qs(self):
        """Always return distinct queryset."""
        parent = super().qs
        return parent.distinct()


class VariantAnnotationFilter(CheckQueryParamsMixin, filters.FilterSet):
    """Filter the VariantAnnotation objects endpoint."""

    class Meta:
        """Filter configuration."""

        model = VariantAnnotation
        fields = {
            "id": NUMBER_LOOKUPS,
            "variant": RELATED_LOOKUPS,
            "type": TEXT_LOOKUPS,
            "transcripts__annotation": TEXT_LOOKUPS,
            "transcripts__annotation_impact": TEXT_LOOKUPS,
            "transcripts__gene": TEXT_LOOKUPS,
            "transcripts__protein_impact": TEXT_LOOKUPS,
            "clinical_diagnosis": TEXT_LOOKUPS,
            "clinical_significance": TEXT_LOOKUPS,
            "dbsnp_id": TEXT_LOOKUPS,
            "clinvar_id": TEXT_LOOKUPS,
        }


class VariantCallFilter(CheckQueryParamsMixin, filters.FilterSet):
    """Filter the VariantCall objects endpoint."""

    class Meta:
        """Filter configuration."""

        model = VariantCall
        fields = {
            "id": NUMBER_LOOKUPS,
            "sample__slug": TEXT_LOOKUPS,
            "sample": RELATED_LOOKUPS,
            "data__slug": TEXT_LOOKUPS,
            "data": RELATED_LOOKUPS,
            "variant": RELATED_LOOKUPS,
            "variant__species": TEXT_LOOKUPS,
            "variant__genome_assembly": TEXT_LOOKUPS,
            "variant__chromosome": TEXT_LOOKUPS,
            "variant__position": NUMBER_LOOKUPS,
            "variant__reference": TEXT_LOOKUPS,
            "variant__alternative": TEXT_LOOKUPS,
            "experiment": RELATED_LOOKUPS,
            "quality": NUMBER_LOOKUPS,
            "depth": NUMBER_LOOKUPS,
            "filter": TEXT_LOOKUPS,
            "genotype": TEXT_LOOKUPS,
            "genotype_quality": NUMBER_LOOKUPS,
            "alternative_allele_depth": NUMBER_LOOKUPS,
            "depth_norm_quality": NUMBER_LOOKUPS,
        }


class VariantExperimentFilter(CheckQueryParamsMixin, filters.FilterSet):
    """Filter the VariantExperiment objects endpoint."""

    class Meta:
        """Filter configuration."""

        model = VariantExperiment
        fields = {
            "id": NUMBER_LOOKUPS,
            "contributor__username": TEXT_LOOKUPS,
            "contributor__email": TEXT_LOOKUPS,
            "contributor__id": NUMBER_LOOKUPS,
            "timestamp": DATETIME_LOOKUPS,
            "variant_data_source": TEXT_LOOKUPS,
        }
