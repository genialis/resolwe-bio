"""
.. autoclass:: server.tests.processor.utils.PreparedData
   :members:

"""
import os
import shutil

from django.conf import settings
from server.models import Data, Processor, iterate_fields
from server.tasks import manager


class PreparedData(object):

    """Mixin to prepare processor test data.

    This class is Mixin can be used with
    :class:`BaseProcessorTestCase` and includes functions for processing
    common data types as inputs for more advanced tests.

    """

    def __init__(self):
        if not hasattr(self, 'case'):
            raise TypeError("PreparedData mixin not used with BaseProcessorTestCase.")

        if not hasattr(self, 'admin'):
            raise TypeError("PreparedData mixin not used with BaseProcessorTestCase.")

        if not hasattr(self, 'assertStatus'):
            raise TypeError("PreparedData mixin not used with BaseProcessorTestCase.")

        if not hasattr(self, 'current_path'):
            raise TypeError("PreparedData mixin not used with BaseProcessorTestCase.")

        self.case = getattr(self, 'case')
        self.admin = getattr(self, 'admin')
        self.assertStatus = getattr(self, 'assertStatus')  # pylint: disable=invalid-name
        self.current_path = getattr(self, 'current_path')

    def run_processor(self, processor_name, input_, assert_status=Data.STATUS_DONE):
        """Runs given processor with specified inputs.

        If input is file, file path should be given relative to
        ``server/tests/processor/inputs`` folder.
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
            author_id=str(self.admin.pk),
            processor_name=processor_name,
            case_ids=[self.case.pk],
        )
        d.save()
        manager(run_sync=True, verbosity=0)

        # Fetch latest Data object from database
        d = Data.objects.get(pk=d.pk)

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
