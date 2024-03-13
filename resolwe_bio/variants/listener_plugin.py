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

        variants = []
        variant_calls = []

        for variant_data in variants_data:
            chromosome = variant_data["chromosome"]
            position = variant_data["position"]
            reference = variant_data["reference"]
            alternative = variant_data["alternative"]
            Variant.objects.exists(
                species=species,
                genome_assembly=genome_assembly,
                chromosome=chromosome,
                position=position,
                reference=reference,
                alternative=alternative,
            )
            variants.append(
                Variant(
                    species=species,
                    genome_assembly=genome_assembly,
                    chromosome=chromosome,
                    position=position,
                    reference=reference,
                    alternative=alternative,
                    type=variant_data.get("type", None),
                )
            )
            quality = variant_data["quality"]
            depth = variant_data["depth"]
            genotype = variant_data["genotype"]
            filter = variant_data["filter"]
            variant_calls.append(
                VariantCall(
                    data=data,
                    sample=sample,
                    quality=quality,
                    depth=depth,
                    genotype=genotype,
                    filter=filter,
                )
            )
        Variant.objects.bulk_create(variants)
