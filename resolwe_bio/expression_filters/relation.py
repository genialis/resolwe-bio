""".. Ignore pydocstyle D400.

=======================
Relations Template Tags
=======================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

from resolwe.flow.expression_engines.jinja.filters import id_
from resolwe.flow.models.entity import RelationPartition

from resolwe_bio.models import Sample


def replicate_groups(data):
    """Get replicate groups."""
    if not isinstance(data, list):
        raise ValueError('List of data IDs must be given')

    data_ids = [id_(d) for d in data]
    samples = Sample.objects.filter(data__id__in=data_ids)

    if len(samples) != len(data_ids):
        raise ValueError('Can not get replicates for data without sample')

    partitions = RelationPartition.objects.filter(relation__category='Replicate', relation__type__name='group',
                                                  entity__in=samples)

    sample_group = {}
    for p in partitions:
        if p.entity.id in sample_group:
            raise ValueError('Can not handle samples with more then one replicate relations')

        sample_group[p.entity.id] = p.label

    group_map = {}
    group_index = 1
    repl_groups = []
    # Ensure the correct order
    for d in data_ids:
        # This is slow because we are fetching samples one by one
        sample_id = Sample.objects.filter(data__id=d).values('id').first()['id']

        if sample_id in sample_group:
            group_label = sample_group[sample_id]

            if group_label not in group_map:
                group_map[group_label] = group_index
                group_index += 1

            repl_groups.append(group_map[group_label])
        else:
            repl_groups.append(group_index)
            group_index += 1

    return repl_groups


# A dictionary of filters that will be registered.
filters = {  # pylint: disable=invalid-name
    'replicate_groups': replicate_groups,
}
