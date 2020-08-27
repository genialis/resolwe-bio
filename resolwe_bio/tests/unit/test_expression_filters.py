from jinja2 import Environment, FunctionLoader

from resolwe.flow.models import Collection, Data, DescriptorSchema, Process, Relation
from resolwe.flow.models.entity import RelationPartition, RelationType
from resolwe.test import TestCase

from resolwe_bio.expression_filters.relation import background_pairs, replicate_groups


class TestReplicateGroupFilter(TestCase):
    def test_repl_groups(self):
        p = Process.objects.create(
            type="data:test:process",
            slug="test-process",
            contributor=self.contributor,
            entity_type="sample",
            entity_descriptor_schema="sample",
        )

        DescriptorSchema.objects.create(
            slug="sample", version="1.0.0", contributor=self.contributor
        )

        collection = Collection.objects.create(
            name="Test collection", contributor=self.contributor
        )

        data = []
        for i in range(10):
            data.append(
                Data.objects.create(
                    name="Data {}".format(i),
                    contributor=self.contributor,
                    process=p,
                    status=Data.STATUS_DONE,
                )
            )

        rel_type_group = RelationType.objects.get(name="group")
        replicate_group = Relation.objects.create(
            contributor=self.contributor,
            collection=collection,
            type=rel_type_group,
            category="Replicate",
        )

        RelationPartition.objects.create(
            relation=replicate_group, entity=data[0].entity, label="Group2"
        )
        RelationPartition.objects.create(
            relation=replicate_group, entity=data[1].entity, label="Group2"
        )
        RelationPartition.objects.create(
            relation=replicate_group, entity=data[2].entity, label="Group2"
        )
        RelationPartition.objects.create(
            relation=replicate_group, entity=data[3].entity, label="Group1"
        )
        RelationPartition.objects.create(
            relation=replicate_group, entity=data[4].entity, label="Group1"
        )
        RelationPartition.objects.create(
            relation=replicate_group, entity=data[6].entity, label="X"
        )
        RelationPartition.objects.create(
            relation=replicate_group, entity=data[7].entity, label="X"
        )
        RelationPartition.objects.create(
            relation=replicate_group, entity=data[8].entity, label="X"
        )

        # Test replicate groups order
        self.assertEqual(
            replicate_groups([{"__id": d.id} for d in data]),
            [1, 1, 1, 2, 2, 3, 4, 4, 4, 5],
        )

        def load_templates(template_name):
            if template_name == "replicate_groups":
                return "{{ some_data|replicate_groups }}"

        env = Environment(loader=FunctionLoader(load_templates))
        env.filters["replicate_groups"] = replicate_groups
        replicate_groups_template = env.get_template("replicate_groups")
        self.assertEqual(
            replicate_groups_template.render(some_data=[{"__id": d.id} for d in data]),
            "[1, 1, 1, 2, 2, 3, 4, 4, 4, 5]",
        )

        # Test list must be given
        with self.assertRaises(ValueError):
            replicate_groups(42)

        collection2 = Collection.objects.create(
            name="Test collection 2", contributor=self.contributor
        )

        replicate_group2 = Relation.objects.create(
            contributor=self.contributor,
            collection=collection2,
            type=rel_type_group,
            category="Replicate",
        )

        RelationPartition.objects.create(
            relation=replicate_group2, entity=data[0].entity, label="Group3"
        )

        # Test ValueError if two relations
        with self.assertRaises(ValueError):
            replicate_groups([{"__id": d.id} for d in data])

        with self.assertRaises(ValueError):
            replicate_groups_template.render(some_data=[{"__id": d.id} for d in data])

        # Test ValueError if repeated data objects
        with self.assertRaises(ValueError):
            replicate_groups([{"__id": data[0].id}, {"__id": data[0].id}])

    def test_repl_groups_no_sample(self):
        p = Process.objects.create(
            type="data:test:process",
            slug="test-process",
            contributor=self.contributor,
        )

        d = Data.objects.create(
            name="Data",
            contributor=self.contributor,
            process=p,
            status=Data.STATUS_DONE,
        )

        with self.assertRaises(ValueError):
            replicate_groups([{"__id": d.id}])


class TestBackgroundPairsFilter(TestCase):
    def test_background_pairs(self):
        proc = Process.objects.create(
            type="data:test:process:",
            slug="test-process",
            contributor=self.contributor,
            entity_type="sample",
            entity_descriptor_schema="sample",
        )

        proc2 = Process.objects.create(
            type="data:test:process2:",
            slug="test-process2",
            contributor=self.contributor,
            entity_type="sample",
            entity_descriptor_schema="sample",
            input_schema=[{"name": "src", "type": "data:test:process:"}],
        )

        DescriptorSchema.objects.create(
            slug="sample", version="1.0.0", contributor=self.contributor
        )

        collection = Collection.objects.create(
            name="Test collection", contributor=self.contributor
        )

        data = []
        data2 = []
        for i in range(6):
            data.append(
                Data.objects.create(
                    name="Data {}".format(i),
                    contributor=self.contributor,
                    process=proc,
                    status=Data.STATUS_DONE,
                )
            )

            data2.append(
                Data.objects.create(
                    name="Data2 {}".format(i),
                    contributor=self.contributor,
                    process=proc2,
                    status=Data.STATUS_DONE,
                    input={"src": data[-1].id},
                )
            )

        rel_type_background = RelationType.objects.get(name="background")
        background = Relation.objects.create(
            contributor=self.contributor,
            collection=collection,
            type=rel_type_background,
            category="Background",
        )

        RelationPartition.objects.create(
            relation=background, entity=data2[0].entity, label="background"
        )
        RelationPartition.objects.create(
            relation=background, entity=data2[1].entity, label="case"
        )
        RelationPartition.objects.create(
            relation=background, entity=data2[3].entity, label="case"
        )
        RelationPartition.objects.create(
            relation=background, entity=data2[4].entity, label="case"
        )

        def load_templates(template_name):
            if template_name == "background_pairs":
                return "{{ some_data|background_pairs }}"

        env = Environment(loader=FunctionLoader(load_templates))
        env.filters["background_pairs"] = background_pairs
        background_pairs_template = env.get_template("background_pairs")

        self.assertEqual(
            background_pairs_template.render(
                some_data=[{"__id": d.id, "__type": d.process.type} for d in data]
            ),
            "[({1}, {0}), ({2}, None), ({3}, {0}), ({4}, {0}), ({5}, None)]".format(
                data[0].id, data[1].id, data[2].id, data[3].id, data[4].id, data[5].id
            ),
        )
        self.assertEqual(
            background_pairs_template.render(
                some_data=[{"__id": d.id, "__type": d.process.type} for d in data2]
            ),
            "[({1}, {0}), ({2}, None), ({3}, {0}), ({4}, {0}), ({5}, None)]".format(
                data2[0].id,
                data2[1].id,
                data2[2].id,
                data2[3].id,
                data2[4].id,
                data2[5].id,
            ),
        )

        # Test list must be given
        with self.assertRaises(ValueError):
            background_pairs(42)
