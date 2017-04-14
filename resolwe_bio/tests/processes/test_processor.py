# pylint: disable=missing-docstring
from django.test import TestCase

from resolwe.flow.models import Process
from resolwe.flow.utils import iterate_schema


class ProcessFieldsTestCase(TestCase):
    def test_processor_types(self):
        procs = list(Process.objects.all())
        types = {}
        errors_equals = set()
        errors_subtype = set()

        for p in procs:
            fields = sorted('{} {}'.format(pth, schema['type']) for schema, _, pth in
                            iterate_schema({}, p.output_schema, 'output'))
            if p.type not in types:
                types[p.type] = {
                    'fields': fields,
                    'name': [p.name]
                }
            else:
                types[p.type]['name'].append(p.name)

                if types[p.type]['fields'] != fields:
                    errors_equals.add(p.type)

        if errors_equals:
            self.fail('Processes of the same type should have the same output fields:\n\n    {}'.format(
                '\n    '.join(', '.join(types[typ]['name']) for typ in errors_equals)))

        type_list = sorted(types)
        for i, typ in enumerate(type_list):
            for prev_typ in type_list[:i]:
                if typ.startswith(prev_typ):
                    prev_typ_fields = types[prev_typ]['fields']
                    typ_fields = types[typ]['fields']
                    if set(prev_typ_fields).difference(typ_fields):
                        errors_subtype.add('{} {}'.format(prev_typ, typ))

        if errors_subtype:
            self.fail('Processors should include all output fields of the parent type:\n\n    {}'.format(
                '\n    '.join(errors_subtype)))
