# Generated by Django 2.2.12 on 2020-07-28 08:06

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("resolwe_bio_kb", "0010_update_indexes_constraints"),
    ]

    operations = [
        migrations.AddIndex(
            model_name="feature",
            index=models.Index(fields=["type"], name="idx_feature_type"),
        ),
    ]
