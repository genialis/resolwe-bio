"""Resolwe-bio process with KnowledgeBase extensions."""

from typing import Any, Dict, List

from resolwe.process.communicator import communicator
from resolwe.process.models import Data
from resolwe.process.runtime import Process

from resolwe_bio.process.models import Feature, Mapping


class ProcessBio(Process):
    """Resolwe-bio process class."""

    # Mark this process abstract ie the base class for real processes.
    _abstract = True

    def __init__(self, data: Data):
        """Initialize."""
        super().__init__(data)
        self.feature = Feature
        self.mapping = Mapping

    def add_variants(self, data_source: str, variants: List[Dict[str, Any]]):
        """Add variants to the database."""
        return communicator.add_variants((data_source, variants))

    def add_variants_annotations(self, variants: List[Dict[str, Any]]):
        """Add variants annotations to the database."""
        return communicator.add_variants_annotations(variants)
