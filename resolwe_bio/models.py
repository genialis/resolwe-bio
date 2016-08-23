""".. Ignore pydocstyle D400.

======
Models
======

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from django.db import models

from resolwe.flow.models import BaseCollection, Collection


class Sample(BaseCollection):
    """Postgres model for storing sample."""

    class Meta(BaseCollection.Meta):
        """Collection Meta options."""

        permissions = (
            ("view_sample", "Can view sample"),
            ("edit_sample", "Can edit sample"),
            ("share_sample", "Can share sample"),
            ("download_sample", "Can download files from sample"),
            ("add_sample", "Can add data objects to sample"),
            ("owner_sample", "Is owner of the sample"),
        )

    #: list of collections to which sample belong
    collections = models.ManyToManyField(Collection)

    #: sample not finalized (missing data, annotations or not yet confirmed)
    presample = models.BooleanField(default=True)
