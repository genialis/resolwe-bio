""".. Ignore pydocstyle D400.

====================
Sample Template Tags
====================

"""
from resolwe_bio.models import Sample


def get_sample(data):
    """Get sample object."""
    if isinstance(data, dict):
        # `Data` object's id is hydrated as `__id` in expression engine
        data = data["__id"]

    sample_qs = Sample.objects.filter(data=data)

    if not sample_qs.exists():
        return None

    return sample_qs.first()


def get_sample_attr(data, attr):
    """Get ``attr`` attribute of sample object."""
    sample = get_sample(data)

    return getattr(sample, attr, None)


def sample_id(data):
    """Return `pk` of `Sample` that given `Data` object belongs to."""
    return get_sample_attr(data, "pk")


def sample_slug(data):
    """Return `slug` of `Sample` that given `Data` object belongs to."""
    return get_sample_attr(data, "slug")


def sample_name(data):
    """Return `name` of `Sample` that given `Data` object belongs to."""
    return get_sample_attr(data, "name")


# A dictionary of filters that will be registered.
filters = {
    "sample_id": sample_id,
    "sample_slug": sample_slug,
    "sample_name": sample_name,
    "sample": get_sample,
}
