from django.db import migrations

from resolwe.flow.migration_ops import ResolweProcessChangeType


class Migration(migrations.Migration):
    """
    Change the ``alignment-star-index`` process type.
    """

    dependencies = [
        ("resolwe_bio", "0012_full_text_search"),
    ]

    operations = [
        ResolweProcessChangeType(
            process="alignment-star-index",
            new_type="data:index:star:",
        ),
    ]
