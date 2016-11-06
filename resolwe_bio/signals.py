""".. Ignore pydocstyle D400.

===============
Signal Handlers
===============

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from django.db.models.signals import pre_delete, post_save
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
        for field_schema, fields, _ in iterate_schema(instance.input, instance.process.input_schema, ''):
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

            for permission in list(zip(*sample._meta.permissions))[0]:  # pylint: disable=protected-access
                shortcuts.assign_perm(permission, sample.contributor, sample)

            # XXX: This doesn't work, because signal is triggered before Data
            #      object is added to collections.
            # for collection in Collection.objects.filter(data=instance.pk):
            #     sample.collections.add(collection)

        sample.data.add(instance)


@receiver(pre_delete, sender=Data)
def add_pre_delete_handler(sender, instance, **kwargs):
    """Delete Sample when last Data object is deleted."""
    try:
        sample = Sample.objects.get(data=instance.pk)
    except Sample.DoesNotExist:  # pylint: disable=no-member
        return

    if sample.data.count() == 1:  # last Data object will be just deleted
        sample.delete()
