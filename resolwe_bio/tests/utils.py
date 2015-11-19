"""
.. autoclass:: server.tests.processor.utils.PreparedData
   :members:

"""
import gzip
import hashlib
import json
import os
import shutil
import zipfile

from django.conf import settings
from django.contrib.auth.models import Group
from django.core import management
from django.test import TestCase

from genesis.models import GenUser
from server.models import Data, Case, Processor, Storage, Template, Trigger, dict_dot, iterate_fields
from server.tasks import manager


PROCESSORS_FIXTURE_CACHE = None


def clear_all():
    """Delete all objects from MongoDB collections."""
    Data.objects.all().delete()
    Case.objects.all().delete()
    Template.objects.all().delete()
    Trigger.objects.all().delete()
    Storage.objects.all().delete()
    GenUser.objects.all().delete()
    Processor.objects.all().delete()


def create_admin(is_superuser=True):
    """Create admin user."""
    clear_all()

    g = Group(name='admins')
    g.save()

    user = GenUser.objects.create_user('admin@admin.admin', 'test_admin_password')
    user.first_name = "Admin"
    user.last_name = "Genialis"
    user.is_superuser = is_superuser
    user.groups.add(g)
    user.save()
    user.password = 'test_admin_password'  # bypass for the hash
    return user


def create_test_case(author_id=None):
    """Create 2 test projects."""
    if not author_id:
        author_id = str(GenUser.objects.first().id)

    ret = {}

    c = Case()
    c.author_id = str(author_id)
    c.name = 'test_case'
    c.description = 'test_description'
    c.url_slug = 'case_slug'
    c.save()
    ret['c1'] = c

    c = Case()
    c.author_id = str(author_id)
    c.name = 'test_case2'
    c.description = 'test_description2'
    c.url_slug = 'new_slug'
    c.save()
    ret['c2'] = c

    return ret


def _register_processors():
    """Register processors.

    Processor definitions are red when the first test is callled and cached
    into the PROCESSORS_FIXTURE_CACHE global variable.

    """
    Processor.objects.delete()

    global PROCESSORS_FIXTURE_CACHE  # pylint: disable=global-statement
    if PROCESSORS_FIXTURE_CACHE:
        Processor.objects.insert(PROCESSORS_FIXTURE_CACHE)
    else:
        if len(GenUser.objects.filter(is_superuser=True)) == 0:
            GenUser.objects.create_superuser(email='admin@genialis.com')

        management.call_command('register', force=True, testing=True, verbosity='0')
        PROCESSORS_FIXTURE_CACHE = Processor.objects.all()
        for p in PROCESSORS_FIXTURE_CACHE:
            # Trick Mongoengine not to fail the insert
            p._created = True  # pylint: disable=protected-access


class ProcessTestCase(TestCase):

    """Base class for writing processor tests.

    This class is subclass of Django's ``TestCase`` with some specific
    functions used for testing processors.

    To write a processor test use standard Django's syntax for writing
    tests and follow next steps:

    #. Put input files (if any) in ``server/tests/processor/files``
       folder.
    #. Run test with run_processor.
    #. Check if processor has finished successfully with
       assertDone function.
    #. Assert processor's output with :func:`assertFiles`,
       :func:`assertFields` and :func:`assertJSON` functions.

    .. DANGER::
        If output files doesn't exists in
        ``tests/files`` folder, they are created
        automatically. But you have to chack that they are correct
        before using them for further runs.

    """

    def setUp(self):
        super(ProcessTestCase, self).setUp()
        self.admin = create_admin()
        _register_processors()

        cases = create_test_case(self.admin.pk)
        self.case = cases['c1']
        self.case2 = cases['c2']
        self.current_path = os.path.dirname(os.path.abspath(__file__))
        self._keep_all = False
        self._keep_failed = False

    def tearDown(self):
        super(ProcessTestCase, self).tearDown()

        # Delete Data objects and their files unless keep_failed
        for d in Data.objects.all():
            if self._keep_all or (self._keep_failed and d.status == "error"):
                print("KEEPING DATA: {}".format(d.pk))
            else:
                data_dir = os.path.join(settings.DATAFS['data_path'], str(d.pk))
                d.delete()
                shutil.rmtree(data_dir, ignore_errors=True)

    def keep_all(self):
        """Do not delete output files after test for all data."""
        self._keep_all = True

    def keep_failed(self):
        """Do not delete output files after test for failed data."""
        self._keep_failed = True

    def assertStatus(self, obj, status):  # pylint: disable=invalid-name
        """Check if Data object's status is 'status'.

        :param obj: Data object for which to check status
        :type obj: :obj:`server.models.Data`
        :param status: Data status to check
        :type status: str

        """
        self.assertEqual(obj.status, status,
                         msg="Data status is '{}', not '{}'".format(obj.status, status) + self._msg_stdout(obj))

    def assertFields(self, obj, path, value):  # pylint: disable=invalid-name
        """Compare Data object's field to given value.

        :param obj: Data object with field to compare
        :type obj: :obj:`server.models.Data`

        :param path: Path to field in Data object.
        :type path: :obj:`str`

        :param value: Desired value.
        :type value: :obj:`str`

        """
        field = dict_dot(obj['output'], path)
        self.assertEqual(field, value,
                         msg="Field 'output.{}' mismatch: {} != {}".format(path, field, str(value)) +
                         self._msg_stdout(obj))

    def assertFiles(self, obj, field_path, fn, compression=None):  # pylint: disable=invalid-name
        """Compare output file of a processor to the given correct file.

        :param obj: Data object which includes file that we want to
            compare.
        :type obj: :obj:`server.models.Data`

        :param field_path: Path to file name in Data object.
        :type field_path: :obj:`str`

        :param fn: File name (and relative path) of file to which we
            want to compare. Name/path is relative to
            'tests/files'.
        :type fn: :obj:`str`

        :param compression: If not None, files will be uncompressed with
            the appropriate compression library before comparison.
            Currently supported compression formats are "gzip" and
            "zip".
        :type compression: :obj:`str`

        """
        if compression is None:
            open_fn = open
        elif compression == 'gzip':
            open_fn = gzip.open
        elif compression == 'zip':
            open_fn = zipfile.ZipFile.open
        else:
            raise ValueError("Unsupported compression format.")

        field = dict_dot(obj['output'], field_path)
        output = os.path.join(settings.DATAFS['data_path'], str(obj.pk), field['file'])
        output_file = open_fn(output)
        output_hash = hashlib.sha256(output_file.read()).hexdigest()

        wanted = os.path.join(self.current_path, 'files', fn)

        if not os.path.isfile(wanted):
            shutil.copyfile(output, wanted)
            self.fail(msg="Output file {} missing so it was created.".format(fn))

        wanted_file = open_fn(wanted)
        wanted_hash = hashlib.sha256(wanted_file.read()).hexdigest()
        self.assertEqual(wanted_hash, output_hash,
                         msg="File hash mismatch: {} != {}".format(wanted_hash, output_hash) + self._msg_stdout(obj))

    def assertFileExist(self, obj, field_path):  # pylint: disable=invalid-name
        """Compare output file of a processor to the given correct file.

        :param obj: Data object which includes file that we want to
            compare.
        :type obj: :obj:`server.models.Data`

        :param field_path: Path to file name in Data object.
        :type field_path: :obj:`str`

        """
        field = dict_dot(obj['output'], field_path)
        output = os.path.join(settings.DATAFS['data_path'], str(obj.pk), field['file'])

        if not os.path.isfile(output):
            self.fail(msg="File {} does not exist.".format(field_path))

    def assertJSON(self, obj, storage, field_path, file_name):  # pylint: disable=invalid-name
        """Compare JSON in Storage object to the given correct output.

        :param obj: Data object which includes file that we want to
            compare.
        :type obj: :obj:`server.models.Data`

        :param storage: Storage (or storage id) which contains JSON to
            compare.
        :type storage: :obj:`server.models.Storage` or :obj:`str`

        :param field_path: Path to JSON subset to compare in Storage
            object. If it is empty, entire Storage object will be
            compared.
        :type field_path: :obj:`str`

        :param file_name: File name (and relative path) of file to which we
            want to compare. Name/path is relative to
            'tests/files'.
        :type file_name: :obj:`str`

        """
        self.assertEqual(os.path.splitext(file_name)[1], '.gz', msg='File extension must be .gz')

        if not isinstance(storage, Storage):
            storage = Storage.objects.get(pk=str(storage))

        storage_obj = dict_dot(storage['json'], field_path)

        file_path = os.path.join(self.current_path, 'files', file_name)
        if not os.path.isfile(file_path):
            with gzip.open(file_path, 'w') as f:
                json.dump(storage_obj, f)

            self.fail(msg="Output file {} missing so it was created.".format(file_name))

        with gzip.open(file_path) as f:
            file_obj = json.load(f)

        self.assertEqual(storage_obj, file_obj,
                         msg="Storage {} field '{}' does not match file {}".format(
                             storage.id, field_path, file_name) + self._msg_stdout(obj))

    def _msg_stdout(self, data):
        """Print stdout.txt content."""
        msg = "\n\nDump stdout.txt:\n\n"
        stdout = os.path.join(settings.DATAFS['data_path'], str(data.pk), 'stdout.txt')
        if os.path.isfile(stdout):
            with open(stdout, 'r') as fn:
                msg += fn.read()

        return msg

    def run_processor(self, processor_name, input_={}, assert_status=Data.STATUS_DONE, run_manager=True, verbosity=0):
        """Runs given processor with specified inputs.

        If input is file, file path should be given relative to
        ``tests/files`` folder.
        If ``assert_status`` is given check if Data object's status
        matches ``assert_status`` after finishing processor.

        :param processor_name: name of the processor to run
        :type processor_name: :obj:`str`

        :param ``input_``: Input paramaters for processor. You don't
            have to specifie parameters for which default values are
            given.
        :type ``input_``: :obj:`dict`

        :param ``assert_status``: Desired status of Data object
        :type ``assert_status``: :obj:`str`

        :return: :obj:`server.models.Data` object which is created by
            the processor.

        """
        p = Processor.objects.get(name=processor_name)

        for field_schema, fields in iterate_fields(input_, p['input_schema']):
            # copy referenced files to upload dir
            if field_schema['type'] == "basic:file:":
                old_path = os.path.join(self.current_path, 'files', fields[field_schema['name']])
                shutil.copy2(old_path, settings.FILE_UPLOAD_TEMP_DIR)
                file_name = os.path.basename(fields[field_schema['name']])
                fields[field_schema['name']] = {
                    'file': file_name,
                    'file_temp': file_name,
                }

            # convert primary keys to strings
            if field_schema['type'].startswith('data:'):
                fields[field_schema['name']] = str(fields[field_schema['name']])
            if field_schema['type'].startswith('list:data:'):
                fields[field_schema['name']] = [str(obj) for obj in fields[field_schema['name']]]

        d = Data(
            input=input_,
            author_id=str(self.admin.pk),
            processor_name=processor_name,
            case_ids=[self.case.pk],
        )
        d.save()
        if run_manager:
            manager(run_sync=True, verbosity=verbosity)

        # Fetch latest Data object from database
        d = Data.objects.get(pk=d.pk)

        if not run_manager and assert_status == Data.STATUS_DONE:
            assert_status = Data.STATUS_RESOLVING

        if assert_status:
            self.assertStatus(d, assert_status)

        return d

    def prepare_genome(self, fn='genome.fasta.gz'):
        """Prepare genome FASTA."""
        inputs = {'src': fn}
        return self.run_processor('import:upload:genome-fasta', inputs)

    def prepare_reads(self, fn='reads.fastq.gz'):
        """Prepare NGS reads FASTQ."""
        inputs = {'src': fn}
        return self.run_processor('import:upload:reads-fastq', inputs)

    def prepare_bam(self, fn='sp_test.bam'):
        """Prepare alignment BAM."""
        inputs = {'src': fn}
        return self.run_processor('import:upload:mapping-bam', inputs)

    def prepare_annotation(self, fn='sp_test.gtf'):
        """Prepare annotation GTF."""
        inputs = {'src': fn}
        return self.run_processor('import:upload:annotation-gtf', inputs)

    def prepare_expression(self, f_rc='00Hr_rc.tab.gz', f_exp='00Hr_tpm.tab.gz', f_type="TPM"):
        """Prepare expression."""
        inputs = {'rc': f_rc, 'exp': f_exp, 'exp_type': f_type}
        return self.run_processor('import:upload:expression', inputs)
