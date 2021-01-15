"""Resolwe-bio process with KnowledgeBase extensions."""

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
