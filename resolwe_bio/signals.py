from __future__ import absolute_import, division, print_function, unicode_literals

from django.db.models.signals import post_save
from django.dispatch import receiver

from guardian import shortcuts

from resolwe.flow.models import Data, DescriptorSchema, iterate_schema
from .models import Sample


@receiver(post_save, sender=Data)
def add_post_save_handler(sender, instance, **kwargs):
    """Add object to flow_collection.

    * Only add `Data object` to `Sample` if process has defined
      `flow_collwection` field.
    * Add object to existing `Sample`, if `input objects` that
      belong to `flow collection` (but not necessary all
      `input objects`), are part of the same `Sample`.
    * If `input objects` belong to different `Samples` or do not belong
      to any `Sample`, create new `Sample`.

    Collect IDs of all `input objects`.

    """
    if kwargs['created'] and instance.process.flow_collection:
        input_objects = []
        for field_schema, fields, path in iterate_schema(instance.input, instance.process.input_schema, ''):
            if 'name' in field_schema and 'type' in field_schema and field_schema['name'] in fields:
                field = fields[field_schema['name']]
                if field_schema['type'].startswith('data:'):
                    input_objects.append(field)
                if field_schema['type'].startswith('list:data:'):
                    input_objects.extend(field)

        sample_query = Sample.objects.filter(data__id__in=input_objects).distinct()

        if sample_query.count() == 1:
            sample = sample_query.first()
        else:
            des_schema = DescriptorSchema.objects.get(slug=instance.process.flow_collection)
            sample = Sample.objects.create(
                contributor=instance.contributor,
                descriptor_schema=des_schema,
                name=instance.name,
            )

            for permission in list(zip(*sample._meta.permissions))[0]:
                shortcuts.assign_perm(permission, sample.contributor, sample)

        sample.data.add(instance)
