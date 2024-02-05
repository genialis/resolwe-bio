""".. Ignore pydocstyle D400.

====================
Sample Template Tags
====================

"""

from resolwe.flow.models.annotations import AnnotationField

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


def sample_annotation(data, annotation_path):
    """Return the annotation value for the given annotation path."""
    annotations = get_sample_attr(data, "annotations")
    if annotations is None:
        return None
    group_name, field_name = AnnotationField.group_field_from_path(annotation_path)
    annotation_value = annotations.filter(
        field__name=field_name, field__group__name=group_name
    ).first()
    return annotation_value.value if annotation_value else None


def sample_slug(data):
    """Return `slug` of `Sample` that given `Data` object belongs to."""
    return get_sample_attr(data, "slug")


def sample_name(data):
    """Return `name` of `Sample` that given `Data` object belongs to."""
    return get_sample_attr(data, "name")


# A dictionary of filters that will be registered.
filters = {
    "sample_annotation": sample_annotation,
    "sample_id": sample_id,
    "sample_slug": sample_slug,
    "sample_name": sample_name,
    "sample": get_sample,
}
