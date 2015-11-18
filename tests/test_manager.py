# pylint: disable=missing-docstring
import os

from django.conf import settings
from django.test import override_settings
from django.core.cache import cache

from tastypie.test import TestApiClient

from .base import BaseProcessorTestCase
from .utils import PreparedData
from ..unit.utils import login_user
from server.models import Data, Trigger
from server.tasks import LXC_IP_KEY_PREFIX, LXC_IP_PREFIX, LXC_IP_SUFFIX_MAX, LXC_IP_SUFFIX_MIN


class RunFunctionTestCase(BaseProcessorTestCase, PreparedData):
    def test_filesize(self):
        reads = self.prepare_reads()
        self.assertFields(reads, 'fastq.size', 715)

    def test_filesize_nonexisting_file(self):
        d = self.run_processor('test:nonexisting', assert_status=Data.STATUS_ERROR)
        self.assertFields(d, 'proc.error', 'Referenced file doesn\'t exists (you/cant/find.me)')

    def test_purge(self):
        d = self.run_processor('test:purge')
        data_path = os.path.join(settings.DATAFS['data_path'], str(d.pk))

        allowed_dirs = [
            os.path.join(data_path, 'foodir'),
            os.path.join(data_path, 'foodir', 'bardir'),
        ]
        allowed_files = [
            os.path.join(data_path, 'stdout.txt'),
            os.path.join(data_path, 'foo'),
            os.path.join(data_path, 'foodir', 'bar'),
            os.path.join(data_path, 'foodir', 'bardir', 'kuku'),
        ]

        # Check if existing files/dirs are allowed
        for path, dirs, files in os.walk(data_path):
            for dir_ in dirs:
                d_path = os.path.join(path, dir_)
                self.assertTrue(d_path in allowed_dirs,
                                msg='Purge should delete {}'.format(d_path))
            for fn in files:
                f_path = os.path.join(path, fn)
                self.assertTrue(f_path in allowed_files,
                                msg='Purge should delete {}'.format(f_path))

        # Check if all allowed dirs exist
        for dir_ in allowed_dirs:
            self.assertTrue(os.path.isdir(dir_), msg='{} is missing'.format(dir_))

        # Check if all allowed files exist
        for fn in allowed_files:
            self.assertTrue(os.path.isfile(fn), msg='{} is missing'.format(fn))

    def test_templatetags(self):
        reads = self.prepare_reads()
        d = self.run_processor('test:templatetags', {'reads': reads.pk})
        self.assertFields(d, 'reads_name', 'reads.fastq.gz (Upload)')
        self.assertFields(d, 'reads_id', str(reads.pk))
        self.assertFields(d, 'reads_type', 'data:reads:fastq:single:')

    def test_invalid_templatetag(self):
        d = self.run_processor('test:invalid-templatetag', assert_status=Data.STATUS_ERROR)
        self.assertFields(
            d, 'proc.error',
            'Could not parse the remainder: \'("data:reads:")\' from \'reads.type.startswith("data:reads:")\'')

    def test_failed_processor(self):
        d = self.run_processor('test:failed', assert_status=Data.STATUS_ERROR)
        self.assertNotEqual(d.date_finish, None)


class ManagerFunctionTestCase(BaseProcessorTestCase, PreparedData):
    def setUp(self):
        super(ManagerFunctionTestCase, self).setUp()
        cache.clear()

    def test_process_after_failure(self):
        self.run_processor('test:failed', run_manager=False)
        self.run_processor('test:invalid-templatetag', run_manager=False)
        self.assertEqual(len(Data.objects.filter(status=Data.STATUS_RESOLVING)), 2)

        self.run_processor('test:success')
        self.assertEqual(len(Data.objects.filter(status=Data.STATUS_RESOLVING)), 0)
        self.assertEqual(len(Data.objects.filter(status=Data.STATUS_DONE)), 1)
        self.assertEqual(len(Data.objects.filter(status=Data.STATUS_ERROR)), 2)

    def test_assigning_lxc_ip_addresses(self):
        # pretend that some IP addresses are already taken
        n_ips_taken = 10
        ip_suffix_max = LXC_IP_SUFFIX_MIN + n_ips_taken
        for ip_suffix in range(LXC_IP_SUFFIX_MIN, ip_suffix_max):
            ip_key = LXC_IP_KEY_PREFIX + LXC_IP_PREFIX + str(ip_suffix)
            cache.set(ip_key, False, timeout=None)
        d = self.run_processor('test:lxc_ip_addr')
        self.assertFields(d, 'ip_address', LXC_IP_PREFIX + str(ip_suffix_max))
        self.assertEqual(
            cache.get(LXC_IP_KEY_PREFIX + LXC_IP_PREFIX + str(ip_suffix_max), "Not set"),
            "Not set")

    def test_handle_all_lxc_ips_taken(self):
        # pretend that all IP addresses are already taken
        for ip_suffix in range(LXC_IP_SUFFIX_MIN, LXC_IP_SUFFIX_MAX):
            ip_key = LXC_IP_KEY_PREFIX + LXC_IP_PREFIX + str(ip_suffix)
            cache.set(ip_key, False, timeout=None)
        d = self.run_processor('test:lxc_ip_addr', assert_status=Data.STATUS_RESOLVING)
        self.assertFields(d, 'proc.info', 'Waiting for a free executor.')


class HelperFunctionsTestCase(BaseProcessorTestCase, PreparedData):
    def test_gensave(self):
        d = self.run_processor('test:helpers:gen_save')
        self.assertFields(d, 'save', 'So we saved the world together for a while.')

    def test_info(self):
        d = self.run_processor('test:helpers:gen_info')
        self.assertFields(d, 'proc.info', 'Dude, I just lied to a samuri.')

    def test_warning(self):
        d = self.run_processor('test:helpers:gen_warning')
        self.assertFields(d, 'proc.warning', 'You might have been in ******, Doc')

    def test_error(self):
        d = self.run_processor('test:helpers:gen_error', assert_status=Data.STATUS_ERROR)
        self.assertFields(d, 'proc.error', 'Maybe you should do it again, just to be sure.')

    def test_save_file(self):
        d = self.run_processor('test:helpers:gen_save_file')
        self.assertFields(d, 'output_file.file', 'test.txt')
        self.assertFields(d, 'output_file.refs', ['ref.txt'])

    def test_checkrc(self):
        d = self.run_processor('test:helpers:gen_checkrc_1', assert_status=Data.STATUS_ERROR)
        self.assertFields(d, 'proc.error', 'Did you just kill that bunny?')
        self.assertFields(d, 'proc.rc', 42)

        d = self.run_processor('test:helpers:gen_checkrc_2')
        self.assertFields(d, 'proc.error', '')
        self.assertFields(d, 'proc.rc', 0)


class RunTriggersTestCase(BaseProcessorTestCase, PreparedData):
    def setUp(self):
        super(RunTriggersTestCase, self).setUp()

        self.api_client = TestApiClient()

        self.d = self.run_processor('test:triggers:create-file')

        t = Trigger()
        t.name = 'test.trigger'
        t.type = 'data:test:triggers:create:'
        t.author_id = str(self.admin.pk)
        t.trigger_input = 'in_file'
        t.input = {'in_file': 'INPUT'}
        t.processor_name = 'test:triggers:check-file'
        t.case_id = str(self.case2.id)
        t.autorun = True
        t.save()

    @override_settings(TESTING=False)
    def test_copy_before_triggers(self):
        login_user(self.api_client, self.admin)

        resp = self.api_client.post(
            '/api/v1/data/{}/copy/'.format(self.d.id),
            format='json',
            data={'case_ids': [self.case2.id]})
        self.assertEqual(resp.status_code, 201)  # pylint: disable=no-member

        all_objs = Data.objects.filter(case_ids__contains=str(self.case2.pk))
        error_objs = Data.objects.filter(case_ids__contains=str(self.case2.pk),
                                         status=Data.STATUS_ERROR)
        # done_objs = Data.objects.filter(case_ids__contains=str(self.case2.pk),
        #                                 status=Data.STATUS_DONE)
        self.assertEqual(len(all_objs), 2)
        self.assertEqual(len(error_objs), 0)
        # TODO: Add this
        # self.assertEqual(len(done_objs), 2)
