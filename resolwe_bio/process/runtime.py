"""Resolwe-bio process with KnowledgeBase extensions."""

from typing import Any, Dict, List

from resolwe.process.communicator import communicator
from resolwe.process.models import Data
from resolwe.process.runtime import Process

from resolwe_bio.process.models import Feature, Mapping, Variant, VariantCall


class ProcessBio(Process):
    """Resolwe-bio process class."""

    # Mark this process abstract ie the base class for real processes.
    _abstract = True

    def __init__(self, data: Data):
        """Initialize."""
        super().__init__(data)
        self.feature = Feature
        self.mapping = Mapping
        self.variant = Variant
        self.variant_call = VariantCall

    def add_variants(self, data_source: str, variants: List[Dict[str, Any]]):
        """Add variants to the database."""
        return communicator.add_variants((data_source, variants))

    def add_variants_annotations(self, variants: List[Dict[str, Any]]):
        """Add variants annotations to the database."""
        BATCH_SIZE = 1000
        for i in range(0, len(variants), BATCH_SIZE):
            communicator.add_variants_annotations(variants[i : i + BATCH_SIZE])
