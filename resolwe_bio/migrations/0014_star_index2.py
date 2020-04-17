from django.db import migrations

from resolwe.flow.utils.iterators import iterate_schema


def migrate_star_downstream_processes(apps, schema_editor):
    """Migrate schemas of processes that use alignment-star as input."""
    Process = apps.get_model("flow", "Process")

    for process in Process.objects.all():
        for schema, _, _ in iterate_schema({}, process.input_schema):
            if schema["type"] == "data:genomeindex:star:":
                schema["type"] = "data:index:star:"
                process.save()


class Migration(migrations.Migration):
    """Migrate schemas of processes that use alignment-star as input."""

    dependencies = [
        ("resolwe_bio", "0013_star_index"),
    ]

    operations = [
        migrations.RunPython(migrate_star_downstream_processes),
    ]
