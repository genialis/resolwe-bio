import gzip
import os
import shutil

from django.conf import settings
from django.db import migrations

from resolwe.flow.migration_ops import (
    DataDefaultOperation,
    ResolweProcessAddField,
    ResolweProcessRenameField,
)
from resolwe.flow.utils import iterate_fields

FASTA_SCHEMA = {
    "name": "fasta",
    "label": "FASTA file",
    "type": "basic:file:",
}


class DefaultUnzipFasta(DataDefaultOperation):
    """Set default value."""

    def prepare(self, data, from_state):
        pass

    def get_default_for(self, data, from_state):
        """Return default for given data object."""
        fastagz = os.path.join(
            settings.FLOW_EXECUTOR["DATA_DIR"],
            data.location.subpath,
            data.output["fastagz"]["file"],
        )
        assert fastagz.endswith(".gz")
        fasta = fastagz[:-3]

        # Decompress.
        with gzip.open(fastagz, "rb") as infile, open(fasta, "wb") as outfile:
            shutil.copyfileobj(infile, outfile)
        # Change owner of this file to the same owner as fastagz:
        fastagz_stat = os.stat(fastagz)
        os.chown(fasta, fastagz_stat.st_uid, fastagz_stat.st_gid)

        # Set modification time of fai file to "now". This is because
        # bedtools requires that modification time of fai is larger than
        # the one of fasta.
        fai = os.path.join(
            settings.FLOW_EXECUTOR["DATA_DIR"],
            data.location.subpath,
            data.output["fai"]["file"],
        )
        os.utime(fai)

        size = os.path.getsize(fasta)
        return {
            "file": os.path.basename(fasta),
            "size": size,
            "total_size": size,
        }


def recompute_data_size(apps, schema_editor):
    """Recompute size of all data objects of process ``upload-fasta-nucl``."""
    Data = apps.get_model("flow", "Data")
    for data in Data.objects.filter(process__slug="upload-fasta-nucl"):
        hydrate_size(data)
        data.save()


def hydrate_size(data):
    """Compute size of all Data object outputs and its cumultative size.

    This is a simplified version of original ``hydrate_size`` function,
    since we need just a subset of it.
    """

    def add_file_size(obj):
        """Add file size to the basic:file field."""
        path = os.path.join(
            settings.FLOW_EXECUTOR["DATA_DIR"], data.location.subpath, obj["file"]
        )

        obj["size"] = os.path.getsize(path)
        obj["total_size"] = obj["size"]

    data_size = 0
    for field_schema, fields in iterate_fields(data.output, data.process.output_schema):
        name = field_schema["name"]
        value = fields[name]
        if "type" in field_schema:
            if field_schema["type"].startswith("basic:file:"):
                add_file_size(value)
                data_size += value.get("total_size", 0)

    data.size = data_size


class Migration(migrations.Migration):
    """
    Make outputs of ``upload-fasta-nucl`` consistent with ``upload-genome``.

    Process ``upload-genome`` stores compressed output in ``fastagz``
    and uncompressed in ``fasta``. Process ``upload-fasta-nucl`` stores
    compressed output in ``fasta`` output field and does not have a
    field with uncompressed output. Therefore ``fasta`` field is first
    renamed to ``fastagz``. Only then ``fasta`` field is added with
    decompressed content.
    """

    dependencies = [
        ("resolwe_bio", "0010_add_relation_types"),
        ("flow", "0028_add_data_location"),
    ]

    operations = [
        ResolweProcessRenameField(
            process="upload-fasta-nucl", field="output.fasta", new_field="fastagz",
        ),
        ResolweProcessAddField(
            process="upload-fasta-nucl",
            field="output.fasta",
            schema=FASTA_SCHEMA,
            default=DefaultUnzipFasta(),
        ),
        migrations.RunPython(recompute_data_size),
    ]
