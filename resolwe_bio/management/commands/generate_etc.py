""".. Ignore pydocstyle D400.

================================
Generate Expression Time Courses
================================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import gzip
import json
import os
import random
import datetime

from django.core.management.base import BaseCommand
from django.conf import settings
from django.utils import timezone
from resolwe.flow.models import Data, Storage
from .utils import get_descriptorschema, get_process, get_superuser


class Command(BaseCommand):
    """Generate ETC objects."""

    help = "Generate ETC objects"

    def add_arguments(self, parser):
        """Command arguments."""
        parser.add_argument('-e', '--n-etc', type=int, default=3,
                            help="Number of ETC objects to generate (default: %(default)s)")
        parser.add_argument('--rseed', action='store_true', help="Use fixed random seed")

    @staticmethod
    def create_etc(gene_ids, path):
        """Generate an expression time course (ETC).

        Besides returning the JSON dump of the generated ETC, store it as a
        gzipped object to the provided path.

        :return: JSON dump of the generated expression time course
        :rtype: str

        """
        times = (0, 4, 8, 12, 16, 20, 24)
        gene_etcs = {}
        with gzip.open(gene_ids, mode='rt') as gene_ids:
            all_genes = [line.strip() for line in gene_ids]
            for gene in all_genes:
                etc = tuple(round(random.gammavariate(1, 100), 2) for _ in range(len(times)))
                gene_etcs[gene] = etc

        json_dump = json.dumps({'etc': {'genes': gene_etcs, 'timePoints': times}}, indent=4, sort_keys=True)
        with gzip.open(os.path.join(path, 'etc.json.gz'), 'wt') as gzip_file:
            gzip_file.write(json_dump)
        return json_dump

    @staticmethod
    def generate_etc_desciptor():
        """Generate the expression time course descriptor."""
        project = [('1.', 'D. discoideum vs. D. purpureum'),
                   ('2.', 'Filter Development vs. cAMP Pulsing; Frequent Sampling'),
                   ('3.', 'GtaC: WT vs. mutants'),
                   ('4.', 'lncRNA transcriptome')]

        projct_number, project_name = random.choice(project)

        annotation = {'projectNumber': projct_number,
                      'project': project_name,
                      'citation': {'name': 'Rosengarten et. al.',
                                   'url': 'http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-015-1491-7'},
                      'treatment': random.choice(['cAMP Pulses', 'Filter Development']),
                      'parental_strain': 'AX4',
                      'growth': 'K. pneumoniae'}

        return annotation

    def create_data(self, reads_name='seq_reads', rseed=None):
        """Generate expression data."""
        # get test data paths
        data_dir = settings.FLOW_EXECUTOR['DATA_DIR']
        test_files_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', '..', 'tests', 'files'))
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

        etc_json_dump = self.create_etc(dicty_genes, os.path.join(data_dir, str(etc.id)))

        json_object = Storage.objects.create(
            json=json.loads(etc_json_dump),
            contributor=get_superuser(),
            data=etc
        )

        etc.output = {
            'etcfile': {'file': 'etc.json.gz'},
            'etc': json_object.id
        }

        etc.status = Data.STATUS_DONE
        etc.save()

        with open(os.path.join(data_dir, str(etc.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload ETC file. Data object was created '
                         'with the generate_etc django-admin command.')

    def handle(self, *args, **options):
        """Command handle."""
        if options['rseed']:
            random.seed(42)
        for _ in range(options['n_etc']):
            self.create_data()
