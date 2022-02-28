from django.db import migrations

from resolwe.flow.migration_ops import ResolweProcessChangeType


class Migration(migrations.Migration):
    """
    Change the ``upload-orange-metadata`` process type.
    """

    dependencies = [
        ("resolwe_bio", "0001_squashed_0015_sample_indices"),
    ]

    operations = [
        ResolweProcessChangeType(
            process="upload-orange-metadata",
            new_type="data:metadata:unique:",
        ),
    ]
