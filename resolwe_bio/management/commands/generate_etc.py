from __future__ import absolute_import, division, print_function, unicode_literals

import csv
import gzip
import json
import os
import random
import datetime

from django.core.management.base import BaseCommand
from django.contrib.auth import get_user_model
from django.conf import settings
from django.utils import timezone
from resolwe.flow.models import Data, DescriptorSchema, Process, Storage
from resolwe_bio.tools.utils import gzopen
from .utils import get_descriptorschema, get_process, get_superuser


class Command(BaseCommand):

    """Generate ETC objects"""

    help = "Generate ETC objects"

    def add_arguments(self, parser):
        parser.add_argument('-e', '--n-etc', type=int, default=3,
                            help="Number of ETC objects to generate (default: %(default)s)")
        parser.add_argument('--rseed', action='store_true', help="Use fixed random seed")


    @staticmethod
    def create_etc(gene_ids, path, file_name):
        genes = {}
        with open(os.path.join(path, file_name), 'w') as f:
            csvwriter = csv.writer(f, delimiter=str('\t'))
            times = [0, 4, 8, 12, 16, 20, 24]
            csvwriter.writerow(['Gene'] + times)
            with gzopen(gene_ids) as gene_ids:
                all_genes = [line.strip() for line in gene_ids]
                for gene in all_genes:
                    expression = [str(random.gammavariate(1, 100)) for i in range(7)]
                    genes[gene] = expression
                    csvwriter.writerow([gene] + expression)

        with open(os.path.join(path, 'etc.json'), 'w') as json_file:
            js_out = '{{"etc":{}}}'.format(json.dumps({'genes': genes, 'timePoints': times}, separators=(',', ':')))
            json_file.write(js_out)
            gzip.open(os.path.join(path, 'etc.json.gz'), 'wb').write(js_out.encode('utf-8'))


    @staticmethod
    def generate_etc_desciptor():
        project = [('1.', 'D. discoideum vs. D. purpureum'),
                   ('2.', 'Filter Development vs. cAMP Pulsing; Frequent Sampling'),
                   ('3.', 'GtaC: WT vs. mutants'),
                   ('4.', 'lncRNA transcriptome')]

        pn, pr = random.choice(project)

        annotation = {'projectNumber': pn,
                      'project': pr,
                      'citation': {'name': 'Rosengarten et. al.',
                                   'url': 'http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1491-7'},
                      'treatment': random.choice(['cAMP Pulses', 'Filter Development']),
                      'parental_strain': 'AX4',
                      'growth': 'K. pneumoniae'}

        return annotation


    def create_data(self, reads_name='seq_reads', rseed=None):
        # get test data paths
        data_dir = settings.FLOW_EXECUTOR['DATA_DIR']
        test_files_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', '..', 'tests', 'processes', 'files'))
        dicty_genes = os.path.join(test_files_path, 'dicty_genes.tab.gz')

        # Create reads data object
        started = timezone.now()
        etc = Data.objects.create(
            slug='etc',
            name="D. discoideum",
            started=started,
            finished=started + datetime.timedelta(minutes=1),
            descriptor_schema=get_descriptorschema('dicty-etc'),
            descriptor=self.generate_etc_desciptor(),
            status=Data.STATUS_PROCESSING,
            process=get_process(slug='upload-etc'),
            contributor=get_superuser(),
            input={'src': {'file': 'etc.tab'}})

        os.mkdir(os.path.join(data_dir, str(etc.id)))

        self.create_etc(dicty_genes, os.path.join(data_dir, str(etc.id)), 'etc.tab')

        json_object = Storage.objects.create(
            json=json.load(open(os.path.join(data_dir, str(etc.id), 'etc.json'))),
            contributor=get_superuser(),
            data=etc)

        etc.output={
                'etcfile': {'file': 'etc.json.gz'},
                'etc': json_object.id
        }

        etc.status = Data.STATUS_DONE
        etc.save()

        with open(os.path.join(data_dir, str(etc.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload ETC file. Data object was created '
                         'with the generate_etc django-admin command.')

    def handle(self, *args, **options):
        if options['rseed']:
            random.seed(42)
        for i in range(options['n_etc']):
            self.create_data()
