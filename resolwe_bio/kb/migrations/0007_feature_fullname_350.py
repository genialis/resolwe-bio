# -*- coding: utf-8 -*-
# Generated by Django 1.11.14 on 2018-07-12 02:15

from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("resolwe_bio_kb", "0006_feature_fullname_300"),
    ]

    operations = [
        migrations.AlterField(
            model_name="feature",
            name="full_name",
            field=models.CharField(blank=True, max_length=350),
        ),
    ]
