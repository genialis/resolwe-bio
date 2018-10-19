""".. Ignore pydocstyle D400.

=======================
Relations Template Tags
=======================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import copy

from jinja2 import contextfilter

from resolwe.flow.expression_engines.jinja.filters import id_
from resolwe.flow.models.entity import RelationPartition, Relation
from resolwe.flow.models.utils import hydrate_input_references, hydrate_input_uploads

from resolwe_bio.models import Sample


def replicate_groups(data):
    """Get replicate groups."""
    if not isinstance(data, list):
        raise ValueError('List of data must be given')

    data_ids = [id_(d) for d in data]

    if len(data_ids) != len(set(data_ids)):
        raise ValueError('Repeated data objects not allowed')

    samples = Sample.objects.filter(data__id__in=data_ids)

    if len(samples) != len(data_ids):
        raise ValueError('Can not get replicates of data without sample')

    partitions = RelationPartition.objects.filter(
        relation__category='Replicate',
        relation__type__name='group',
        entity__in=samples
    )

    sample_group = {}
    for p in partitions:
        if p.entity.id in sample_group:
            raise ValueError('More than one replicate relation on sample: {}'.format(p.entity.id))

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


@contextfilter
def background_data(context, data_obj):
    """Get background data object"""
    current_data_id = context['proc']['data_id']

    partition = RelationPartition.objects.filter(
        relation__type__name='background',
        label = 'background',
        relation__relationpartition__label='case',
        relation__relationpartition__entity__data=data_obj['__id'],
    )[:2]

    if len(partition) == 0:
        raise ValueError('Background sample not found for data {}'.format(data_obj['__id']))

    if len(partition) > 1:
        raise ValueError('More than one background sample defined for data: {}'.format(data_obj['__id']))

    sample = partition[0].entity
    data = sample.data.filter(process__type=data_obj['__type'])[:2]

    if len(data) == 0:
        raise ValueError('Background sample found but no matching data for data {}'.format(data_obj['__id']))

    if len(data) > 1:
        raise ValueError('Background sample found but 2 or more matching data for data {}'.format(data_obj['__id']))

    data = data[0]
    inputs = copy.deepcopy(data.input)
    hydrate_input_references(inputs, data.process.input_schema)
    hydrate_input_uploads(inputs, data.process.input_schema)

    return inputs


# A dictionary of filters that will be registered.
filters = {  # pylint: disable=invalid-name
    'background_data': background_data,
    'replicate_groups': replicate_groups,
}
