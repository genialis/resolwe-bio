"""Handle variants related commands."""

import logging
import re
from collections import defaultdict
from contextlib import suppress
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Dict, List, Union

from django.conf import settings
from django.db import transaction
from django.db.models import Sum
from django.db.models.functions import Coalesce
from django.urls import reverse
from django.utils.timezone import now

from resolwe.flow.executors.socket_utils import Message, Response
from resolwe.flow.managers.protocol import ExecutorProtocol
from resolwe.flow.models import Data, DataDependency, Process, Worker
from resolwe.flow.models.utils import validate_data_object
from resolwe.flow.utils import dict_dot, iterate_fields, iterate_schema
from resolwe.storage.connectors import connectors
from resolwe.storage.connectors.hasher import StreamHasher
from resolwe.storage.models import ReferencedPath, StorageLocation
from resolwe.utils import BraceMessage as __

from .models import Variant, VariantCall
from .plugin import ListenerPlugin, listener_plugin_manager

if TYPE_CHECKING:
    from resolwe.flow.managers.listener.listener import Processor

logger = logging.getLogger(__name__)


class BasicCommands(ListenerPlugin):
    """Basic listener handlers."""

    plugin_manager = listener_plugin_manager

    def add_variants(
        self, data_id: int, message: Message[dict], manager: "Processor"
    ) -> Response[int]:
        """Handle connecting variants with the samples.

        If the reported variant does not exist in the file it is created.
        """
        data = manager.data(data_id)
        sample = data.sample
        metadata, variants_data = message.message_data
        species, genome_assembly = metadata["species"], metadata["genome_assembly"]

        variant_calls = list()
        variant_cache = dict()

        # Bulk create variants. The consequesce of ignore_conflicts flag is that the
        # database does not returt the ids of the created objects. So first create all
        # the variants and then create the variant calls.
        for variant_data in variants_data:
            key = {
                "species": species,
                "genome_assembly": genome_assembly,
                "chromosome": variant_data["chromosome"],
                "position": variant_data["position"],
                "reference": variant_data["reference"],
                "alternative": variant_data["alternative"],
            }
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
                    quality=variant_data["quality"],
                    depth=variant_data["depth"],
                    genotype=variant_data["genotype"],
                    filter=variant_data["filter"],
                    variant=variant,
                )
            )

        VariantCall.objects.bulk_create(variant_calls)
