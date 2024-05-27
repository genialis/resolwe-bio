"""Handle variants related commands."""

import logging
from typing import TYPE_CHECKING, NotRequired, TypedDict

from resolwe.flow.executors.socket_utils import Message, Response
from resolwe.flow.managers.listener.plugin import (
    ListenerPlugin,
    listener_plugin_manager,
)
from resolwe.flow.models import Data

from .models import (
    Variant,
    VariantAnnotation,
    VariantAnnotationTranscript,
    VariantCall,
    VariantExperiment,
    VariantPrimaryKey,
)

if TYPE_CHECKING:
    from resolwe.flow.managers.listener.listener import Processor

logger = logging.getLogger(__name__)


class VariantAnnotationTranscriptData(TypedDict):
    """Variant annotation transcript data dictionary."""

    annotation: str
    annotation_impact: str
    gene: str
    protein_impact: str
    transcript_ids: list[str]
    canonical: NotRequired[bool]


class VariantAnnotationData(TypedDict):
    """Variant annotation data dictionary."""

    # Variant primary key
    species: str
    genome_assembly: str
    chromosome: str
    position: int
    reference: str
    alternative: str

    # VariantAnnotation data
    type: NotRequired[str]
    clinical_diagnosis: NotRequired[str]
    clinical_significance: NotRequired[str]
    dbsnp_id: NotRequired[str]
    clinvar_id: NotRequired[str]
    data: NotRequired[Data]

    # Transcript data
    transcripts: list[VariantAnnotationTranscriptData]


class VariantData(TypedDict):
    """Variant data dictionary."""

    species: str
    genome_assembly: str
    chromosome: str
    position: int
    reference: str
    alternative: str

    quality: NotRequired[float]
    depth_norm_quality: NotRequired[float]
    unfiltered_allele_depth: NotRequired[int]
    depth: NotRequired[int]
    genotype: NotRequired[str]
    genotype_quality: NotRequired[int]
    filter: NotRequired[str]


class VariantCommands(ListenerPlugin):
    """Listener handlers related to the variants application."""

    plugin_manager = listener_plugin_manager

    def handle_add_variants(
        self,
        data_id: int,
        message: Message[tuple[str, list[VariantData]]],
        manager: "Processor",
    ) -> Response[int]:
        """Handle connecting variants with the samples.

        If the reported variant does not exist in the file it is created.
        """
        # The data object must be included in the sample.
        data = manager.data(data_id)
        sample = data.entity

        if not sample:
            raise RuntimeError("Data object must be included in the sample.")

        # Create an experiment.
        data_source: str = message.message_data[0]
        variants_data: list[VariantData] = message.message_data[1]

        experiment = VariantExperiment.objects.create(
            variant_data_source=data_source,
            contributor=data.contributor,
        )

        # Create the variants.
        variant_calls = list()
        variant_cache = dict()

        # Bulk create variants. The consequesce of ignore_conflicts flag is that the
        # database does not returt the ids of the created objects. So first create all
        # the variants and then create the variant calls.
        for variant_data in variants_data:
            key = {key: variant_data[key] for key in VariantPrimaryKey}
            # To reduce the hits to the database use cache for variants.
            key_tuple = tuple(key.values())
            if key_tuple not in variant_cache:
                variant_cache[key_tuple] = Variant.objects.get_or_create(**key)[0]
            variant = variant_cache[key_tuple]

            variant_calls.append(
                VariantCall(
                    variant=variant,
                    data=data,
                    sample=sample,
                    quality=variant_data.get("quality"),
                    depth_norm_quality=variant_data.get("depth_norm_quality"),
                    unfiltered_allele_depth=variant_data.get("unfiltered_allele_depth"),
                    depth=variant_data.get("depth"),
                    genotype=variant_data.get("genotype"),
                    genotype_quality=variant_data.get("genotype_quality"),
                    filter=variant_data.get("filter"),
                    experiment=experiment,
                )
            )

        VariantCall.objects.bulk_create(variant_calls)
        return message.respond(experiment.id)

    def handle_add_variants_annotations(
        self,
        data_id: int,
        message: Message[list[VariantAnnotationData]],
        manager: "Processor",
    ) -> Response[int]:
        """Add variants annotations."""
        data = manager.data(data_id)

        typed_data: list[VariantAnnotationData] = message.message_data
        # Create the variants.
        variant_cache: dict[tuple, Variant] = dict()

        for annotation_data in typed_data:

            key = {key: annotation_data[key] for key in VariantPrimaryKey}
            # To reduce the hits to the database use cache for variants.
            key_tuple = tuple(key.values())
            if key_tuple not in variant_cache:
                variant_cache[key_tuple] = Variant.objects.get_or_create(**key)[0]
            variant = variant_cache[key_tuple]
            annotation_key = {
                "variant": variant,
                "type": annotation_data.get("type"),
                "clinical_diagnosis": annotation_data.get("clinical_diagnosis"),
                "clinical_significance": annotation_data.get("clinical_significance"),
                "dbsnp_id": annotation_data.get("dbsnp_id"),
                "clinvar_id": annotation_data.get("clinvar_id"),
                "data": data,
            }

            if VariantAnnotation.objects.filter(**annotation_key).exists():
                VariantAnnotation.objects.get(**annotation_key).delete()

            annotation = VariantAnnotation.objects.create(**annotation_key)
            VariantAnnotationTranscript.objects.bulk_create(
                VariantAnnotationTranscript(
                    variant_annotation=annotation,
                    annotation=transcript_data["annotation"],
                    annotation_impact=transcript_data["annotation_impact"],
                    gene=transcript_data["gene"],
                    protein_impact=transcript_data["protein_impact"],
                    transcript_ids=transcript_data["transcript_ids"],
                    canonical=transcript_data.get("canonical", False),
                )
                for transcript_data in annotation_data["transcripts"]
            )
