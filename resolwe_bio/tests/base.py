"""
:class:`server.tests.processor.base.BaseProcessorTestCase`

 .. autoclass:: server.tests.processor.base.BaseProcessorTestCase
    :members: assertDone, assertFields, assertFiles, assertJSON
"""

import hashlib
import gzip
import os
import shutil

from django.conf import settings
from django.core import management
from django.test import TestCase

from server.models import Data, Storage, Processor, GenUser, iterate_fields
from server.tasks import manager
from ..unit.utils import create_test_case, clear_all


class BaseProcessorTestCase(TestCase):
    @classmethod
    def setUpClass(cls):
        clear_all()

        super(BaseProcessorTestCase, cls).setUpClass()

        cls.user = GenUser.objects.create_superuser(email="admin@genialis.com")
        management.call_command('register')

        cls.created_files = []

    @classmethod
    def tearDownClass(cls):
        super(BaseProcessorTestCase, cls).tearDownClass()

        if len(cls.created_files):
            print "#" * 80
            print "WARNING: Next files were created: {}".format(', '.join(cls.created_files))
            print "#" * 80

    def setUp(self):
        super(BaseProcessorTestCase, self).setUp()

        self.current_path = os.path.dirname(os.path.abspath(__file__))
        self.case = create_test_case(self.user.pk)['c1']

    def tearDown(self):
        super(BaseProcessorTestCase, self).tearDown()

        # Delete Data objects and their files
        for d in Data.objects.all():
            data_dir = os.path.join(settings.DATAFS['data_path'], str(d.pk))
            shutil.rmtree(data_dir, ignore_errors=True)
            d.delete()

    def _get_field(self, obj, path):
        if len(path):
            for p in path.split('.'):
                obj = obj[p]
        return obj

    def run_processor(self, processor_name, input_):
        p = Processor.objects.get(name=processor_name)

        for field_schema, fields in iterate_fields(input_, p['input_schema']):
            # copy referenced files to upload dir
            if field_schema['type'] == "basic:file:":
                old_path = os.path.join(self.current_path, 'inputs', fields[field_schema['name']])
                shutil.copy2(old_path, settings.FILE_UPLOAD_TEMP_DIR)
                file_name = os.path.basename(fields[field_schema['name']])
                fields[field_schema['name']] = {
                    'file': file_name,
                    'file_temp': file_name,
                    'is_remote': False,
                }

            # convert primary keys to strings
            if field_schema['type'].startswith('data:'):
                fields[field_schema['name']] = str(fields[field_schema['name']])
            if field_schema['type'].startswith('list:data:'):
                fields[field_schema['name']] = [str(obj) for obj in fields[field_schema['name']]]

        d = Data(
            input=input_,
            author_id=self.user.pk,
            processor_name=processor_name,
            case_ids=[self.case.pk],
        )
        d.save()
        manager(run_sync=True)

        return Data.objects.get(pk=d.pk)

    def assertDone(self, obj):  # pylint: disable=invalid-name
        """Check if Data object's status id 'done'.

        :param obj: Data object for which to check status
        :type obj: :obj:`server.models.Data`

        """
        self.assertEqual(obj.status, 'done')

    def assertFields(self, obj, path, value):  # pylint: disable=invalid-name
        """Compare Data object's field to given value.

        :param obj: Data object with field to compare
        :type obj: :obj:`server.models.Data`

        :param path: Path to field in Data object.
        :type path: :obj:`str`

        :param value: Desired value.
        :type value: :obj:`str`

        """
        field = self._get_field(obj['output'], path)
        return self.assertEqual(field, str(value))

    def assertFiles(self, obj, field_path, fn, gzipped=False):  # pylint: disable=invalid-name
        """Compare output file of a processor to the given correct file.

        :param obj: Data object which includes file that we want to
            compare.
        :type obj: :obj:`server.models.Data`

        :param field_path: Path to file name in Data object.
        :type field_path: :obj:`str`

        :param fn: File name (and relative path) of file to which we
            want to compare. Name/path is relative to
            'server/tests/processor/outputs'.
        :type fn: :obj:`str`

        :param gzipped: If true, file is unziped before comparison.
        :type gzipped: :obj:`bool`

        """
        field = self._get_field(obj['output'], field_path)
        output = os.path.join(settings.DATAFS['data_path'], str(obj.pk), field['file'])
        output_file = gzip.open(output, 'rb') if gzipped else open(output)
        output_hash = hashlib.sha256(output_file.read()).hexdigest()

        wanted = os.path.join(self.current_path, 'outputs', fn)
        if os.path.isfile(wanted):
            wanted_file = gzip.open(wanted, 'rb') if gzipped else open(wanted)
            wanted_hash = hashlib.sha256(wanted_file.read()).hexdigest()
            return self.assertEqual(wanted_hash, output_hash)

        shutil.copyfile(output, wanted)
        self.created_files.append(fn)

    def assertJSON(self, storage, field_path, fn):  # pylint: disable=invalid-name
        """Compare JSON in Storage object to the given correct output.

        :param storage: Storage (or storage id) which contains JSON to
            compare.
        :type storage: :obj:`server.models.Storage` or :obj:`str`

        :param field_path: Path to JSON subset to compare in Storage
            object. If it is empty, entire Storage object will be
            compared.
        :type field_path: :obj:`str`

        :param fn: File name (and relative path) of file to which we
            want to compare. Name/path is relative to
            'server/tests/processor/outputs'.
        :type fn: :obj:`str`

        """
        if type(storage) is not Storage:
            storage = Storage.objects.get(pk=str(storage))

        field = str(self._get_field(storage['json'], field_path))
        field_hash = hashlib.sha256(field).hexdigest()

        wanted = os.path.join(self.current_path, 'outputs', fn)

        if os.path.isfile(wanted):
            wanted_hash = hashlib.sha256(open(wanted).read()).hexdigest()
            return self.assertEqual(wanted_hash, field_hash)

        with open(wanted, 'w') as file_:
            file_.write(field)
        self.created_files.append(fn)
