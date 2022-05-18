from django.db import migrations

from resolwe.flow.migration_ops import ResolweProcessChangeType


class Migration(migrations.Migration):
    """
    Change the ``merge-fastq-single`` and ``merge-fastq-paired`` process type.
    """

    dependencies = [
        ("resolwe_bio", "0003_sample_indices"),
    ]

    operations = [
        ResolweProcessChangeType(
            process="merge-fastq-single",
            new_type="data:mergereads:single:",
        ),
        ResolweProcessChangeType(
            process="merge-fastq-paired",
            new_type="data:mergereads:paired:",
        ),
    ]
