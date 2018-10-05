# pylint: disable=missing-docstring
from jinja2 import Environment, FunctionLoader

from resolwe.flow.models import Data, DescriptorSchema, Collection, Process, Relation
from resolwe.flow.models.entity import RelationPartition, RelationType
from resolwe.test import TestCase

from resolwe_bio.expression_filters.relation import replicate_groups


class TestRelationFilters(TestCase):

    def test_repl_groups(self):
        p = Process.objects.create(
            type='data:test:process',
            slug='test-process',
            contributor=self.contributor,
            flow_collection='sample',
        )

        DescriptorSchema.objects.create(
            slug='sample',
            version='1.0.0',
            contributor=self.contributor
        )

        collection = Collection.objects.create(
            name='Test collection',
            contributor=self.contributor
        )

        data = []
        for i in range(10):
            data.append(Data.objects.create(
                name='Data {}'.format(i),
                contributor=self.contributor,
                process=p,
                status=Data.STATUS_DONE,
            ))

        rel_type_group = RelationType.objects.get(name='group')
        replicate_group = Relation.objects.create(
            contributor=self.contributor,
            collection=collection,
            type=rel_type_group,
            category='Replicate'
        )

        RelationPartition.objects.create(relation=replicate_group, entity=data[0].entity_set.first(), label='Group2')
        RelationPartition.objects.create(relation=replicate_group, entity=data[1].entity_set.first(), label='Group2')
        RelationPartition.objects.create(relation=replicate_group, entity=data[2].entity_set.first(), label='Group2')
        RelationPartition.objects.create(relation=replicate_group, entity=data[3].entity_set.first(), label='Group1')
        RelationPartition.objects.create(relation=replicate_group, entity=data[4].entity_set.first(), label='Group1')
        RelationPartition.objects.create(relation=replicate_group, entity=data[6].entity_set.first(), label='X')
        RelationPartition.objects.create(relation=replicate_group, entity=data[7].entity_set.first(), label='X')
        RelationPartition.objects.create(relation=replicate_group, entity=data[8].entity_set.first(), label='X')

        # Test replicate groups order
        self.assertEqual(replicate_groups([{'__id': d.id} for d in data]), [1, 1, 1, 2, 2, 3, 4, 4, 4, 5])

        def load_templates(template_name):
            if template_name == 'replicate_groups':
                return '{{ some_data|replicate_groups }}'

        env = Environment(loader=FunctionLoader(load_templates))
        env.filters['replicate_groups'] = replicate_groups
        replicate_groups_template = env.get_template('replicate_groups')
        self.assertEqual(replicate_groups_template.render(some_data=[{'__id': d.id} for d in data]),
                         '[1, 1, 1, 2, 2, 3, 4, 4, 4, 5]')

        # Test list must be given
        with self.assertRaises(ValueError):
            replicate_groups(42)

        collection2 = Collection.objects.create(
            name='Test collection 2',
            contributor=self.contributor
        )

        replicate_group2 = Relation.objects.create(
            contributor=self.contributor,
            collection=collection2,
            type=rel_type_group,
            category='Replicate'
        )

        RelationPartition.objects.create(relation=replicate_group2, entity=data[0].entity_set.first(), label='Group3')

        # Test ValueError if two relations
        with self.assertRaises(ValueError):
            replicate_groups([{'__id': d.id} for d in data])

        with self.assertRaises(ValueError):
            replicate_groups_template.render(some_data=[{'__id': d.id} for d in data])

        # Test ValueError if repeated data objects
        with self.assertRaises(ValueError):
            replicate_groups([{'__id': data[0].id}, {'__id': data[0].id}])

    def test_repl_groups_no_sample(self):
        p = Process.objects.create(
            type='data:test:process',
            slug='test-process',
            contributor=self.contributor,
        )

        d = Data.objects.create(
            name='Data',
            contributor=self.contributor,
            process=p,
            status=Data.STATUS_DONE,
        )

        with self.assertRaises(ValueError):
            replicate_groups([{'__id': d.id}])
