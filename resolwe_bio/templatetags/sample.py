from __future__ import absolute_import, division, print_function, unicode_literals

from django import template

from resolwe_bio.models import Sample


register = template.Library()  # pylint: disable=invalid-name


def get_sample_attr(data, attr):
    if isinstance(data, dict):
        # `Data` object's id is hydrated as `__id` in expression engine
        data = data['__id']

    sample_qs = Sample.objects.filter(data=data)

    if not sample_qs.exists():
        return None

    return getattr(sample_qs.first(), attr, None)


@register.filter
def sample_id(data):
    """Return `pk` of `Sample` that given `Data` object belongs to."""

    return get_sample_attr(data, 'pk')


@register.filter
def sample_slug(data):
    """Return `slug` of `Sample` that given `Data` object belongs to."""

    return get_sample_attr(data, 'slug')


@register.filter
def sample_name(data):
    """Return `name` of `Sample` that given `Data` object belongs to."""

    return get_sample_attr(data, 'name')
