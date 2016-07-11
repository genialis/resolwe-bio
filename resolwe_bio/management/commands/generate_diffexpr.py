from __future__ import absolute_import, division, print_function, unicode_literals

import gzip
import json
import logging
import os
import random
import shutil
import string
import datetime
import json
import csv

from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from django.conf import settings
from django.utils import timezone

from resolwe.flow.models import Data, Process, Storage
from resolwe_bio.models import Sample
from resolwe_bio.tools.utils import gzopen
from .utils import get_descriptorschema, get_process, get_superuser, generate_sample_desciptor


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


class Command(BaseCommand):

    """Generate test differential expression data"""

    help = "Generate differential expressions"

    def add_arguments(self, parser):
        parser.add_argument('-n', '--n-diffexps', type=int, default=4,
                            help="Number of differential expressions to generate (default: %(default)s)")
        parser.add_argument('-g', '--group-size', type=int, default=5,
                            help="Number of samples in case/control group per DE (default: %(default)s)")

    @staticmethod
    def get_random_word(length):
        return ''.join(random.choice(string.ascii_lowercase) for _ in range(length))

    def get_name(self, de_name):
        return 'DE_{}_{}'.format(self.get_random_word(4), de_name)

    def create_expressions(self, n):
        expressions = []
        sample_name = self.get_random_word(4)

        for i in range(n):
            cuffquant_file = 'cuffquant_{}.cxb'.format(random.choice([1, 2]))

            # Create expressios
            exp = Data.objects.create(
                name='Smpl_Ex_{}_rep{}'.format(sample_name, i+1),
                process=get_process('upload-cxb'),
                contributor=get_superuser(),
                status=Data.STATUS_PROCESSING,
                input={'src': {'file': cuffquant_file}})

            os.mkdir(os.path.join(self.data_dir, str(exp.id)))
            shutil.copy(os.path.join(self.test_files_path, cuffquant_file), os.path.join(self.data_dir, str(exp.id)))

            exp.output = {
                'cxb': {'file': cuffquant_file},
            }
            exp.status = Data.STATUS_DONE
            exp.save()

            sample = Sample.objects.filter(data=exp)[0]
            sample.presample = False
            sample.descriptor = generate_sample_desciptor('Hs_')
            sample.save()

            with open(os.path.join(self.data_dir, str(exp.id), 'stdout.txt'), 'w') as stdout:
                stdout.write('Upload gene expressions. Sample was created '
                             'with the generate_de django-admin command.')

            logger.info('Created sample: {} (id={})'.format(sample.name, sample.id))
            logger.info('\tData object: (id={})'.format(exp.id))
            expressions.append(exp)

        return expressions

    def create_genome_annotation(self, filename):
        # Create genome annotation
        ann = Data.objects.create(
            name='Annotation_{}'.format(filename.split('.')[0]),
            process=get_process('upload-gtf'),
            contributor=get_superuser(),
            status=Data.STATUS_PROCESSING,
            input={'src': {'file': filename}})

        os.mkdir(os.path.join(self.data_dir, str(ann.id)))

        with gzip.open(os.path.join(self.test_files_path, filename), 'rb') as gzfile:
            with open(os.path.join(self.data_dir, str(ann.id), filename[:-3:]), 'w') as outfile:
                outfile.write(gzfile.read().decode('utf-8'))

        ann.output = {
            'gtf': {'file': filename[:-3:]}
        }
        ann.status = Data.STATUS_DONE
        ann.save()

        with open(os.path.join(self.data_dir, str(ann.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload genome annotation with the generate_De django-admin command.')

        logger.info('Genome annotation created: {} (id={})'.format(filename, ann.id))

        return ann

    def generate_diffexp_data(self, group_size):

        de_name = 'cuffdiff'
        self.data_dir = settings.FLOW_EXECUTOR['DATA_DIR']
        self.test_files_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', '..', 'tests', 'processes', 'files'))
        de_file = os.path.join(self.test_files_path, de_name + '.tab.gz')

        logger.info('---- case samples ----')
        case_samples = self.create_expressions(group_size)

        logger.info('---- control samples ----')
        control_samples = self.create_expressions(group_size)

        logger.info('---- upload annotation ----')

        case_input = [sample.id for sample in case_samples]
        control_input = [sample.id for sample in control_samples]
        genome_annotation_input = self.create_genome_annotation('hg19_chr20_small.gtf.gz')

        de_descriptor = {
            'fc_threshold': 2,
            'fdr_threshold': 0.01,
            'case_label': 'Case group',
            'control_label': 'Control group'
        }

        de_inputs = {
            'case': case_input,
            'control': control_input,
            'annotation': genome_annotation_input.id
        }

        # Create DE data object
        started = timezone.now()

        de_obj = Data.objects.create(
            name=self.get_name(de_name),
            started=started,
            finished=started + datetime.timedelta(minutes=20),
            status=Data.STATUS_PROCESSING,
            descriptor_schema=get_descriptorschema('diff-exp'),
            descriptor=de_descriptor,
            process=get_process(de_name),
            contributor=get_superuser(),
            input=de_inputs)

        # Create data directory and copy files into it
        os.mkdir(os.path.join(self.data_dir, str(de_obj.id)))
        shutil.copy(de_file, os.path.join(self.data_dir, str(de_obj.id)))

        # Get JSON for volcano plot
        with gzip.open(os.path.join(self.test_files_path, '{}.json.gz'.format(de_name)), 'r') as f:
            de_data = json.loads(f.read().decode('utf-8'))

        # TODO: reference on existing true files
        de_obj.output = {
            'diffexp': {'file': de_file},
            'de_data': de_data,
            'diffexp': {'file': de_file},
            'transcript_diff_exp': {'file': de_file},
            'cds_diff_exp': {'file': de_file},
            'tss_group_diff_exp': {'file': de_file},
            'cuffdiff_output': {'file': de_file}
        }

        de_obj.status = Data.STATUS_DONE
        de_obj.save()

        logger.info('---- new differential expression ----')
        logger.info('DE created with id: {}'.format(de_obj.id))

        # Create stdout file
        with open(os.path.join(self.data_dir, str(de_obj.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Differential expression was created with the generate_de django-admin command.')

    def handle(self, *args, **options):
        for i in range(options['n_diffexps']):
            self.generate_diffexp_data(options['group_size'])
