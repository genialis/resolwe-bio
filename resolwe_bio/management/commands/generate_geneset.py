""".. Ignore pydocstyle D400.

=================
Generate Gene Set
=================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import csv
import gzip
import json
import logging
import os
import random
import string
import datetime

from django.core.management.base import BaseCommand
from django.conf import settings
from django.utils import timezone

from resolwe.flow.models import Data, Storage
from resolwe.utils import BraceMessage as __
from .utils import get_descriptorschema, get_process, get_superuser


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


class Command(BaseCommand):
    """Generate gene set data."""

    help = "Generate test gene sets."

    def __init__(self, *args, **kwargs):
        """Set command defaults."""
        super(Command, self).__init__(*args, **kwargs)

        self.data_dir = settings.FLOW_EXECUTOR['DATA_DIR']
        self.test_files_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', '..', 'tests', 'files'))

    def add_arguments(self, parser):
        """Define command arguments."""
        parser.add_argument('-n', '--n-geneset', type=int, default=5,
                            help="Number of gene sets to generate (default: %(default)s)")
        parser.add_argument('--rseed', action='store_true', help="Use fixed random seed")

    @staticmethod
    def get_random_word(length):
        """Generate a random word."""
        return ''.join(random.choice(string.ascii_lowercase) for _ in range(length))

    @staticmethod
    def generate_geneset_file(gene_ids, num, path):
        """Generate gene set file."""
        with gzip.open(os.path.join(path, 'geneset.tab.gz'), 'wt') as f:
            csvwriter = csv.writer(f, delimiter=str('\t'), lineterminator='\n')
            with gzip.open(gene_ids, mode='rt') as gene_ids:
                all_genes = [line.strip() for line in gene_ids]
                geneset = sorted(set([random.choice(all_genes) for _ in range(num)]))
                for gene in geneset:
                    csvwriter.writerow([gene])

        json_dump = json.dumps({'genes': geneset}, indent=4, sort_keys=True)
        with open(os.path.join(path, 'geneset.json'), 'w') as json_file:
            json_file.write(json_dump)

    def create_geneset(self):
        """Create gene set object."""
        started = timezone.now()
        geneset = Data.objects.create(
            name='GeneSet_{}_{}'.format(random.choice(['Hs', 'Mm']), self.get_random_word(3)),
            process=get_process('upload-geneset'),
            contributor=get_superuser(),
            started=started,
            finished=started + datetime.timedelta(seconds=5),
            descriptor_schema=get_descriptorschema('geneset'),
            descriptor={'description': 'Gene set description.'},
            status=Data.STATUS_PROCESSING,
            input={'src': {'file': 'geneset.tab.gz'}, 'source': 'UCSC'})

        mouse_genes = os.path.join(self.test_files_path, 'mouse_genes.tab.gz')

        os.mkdir(os.path.join(self.data_dir, str(geneset.id)))
        self.generate_geneset_file(mouse_genes,
                                   random.randint(15, 150),
                                   os.path.join(self.data_dir, str(geneset.id)))

        json_object = Storage.objects.create(
            json=json.load(open(os.path.join(self.data_dir, str(geneset.id), 'geneset.json'))),
            contributor=get_superuser(),
            name='{}_storage'.format(geneset.name),
            data=geneset)

        os.remove(os.path.join(self.data_dir, str(geneset.id), 'geneset.json'))

        geneset.output = {
            'geneset': {'file': 'geneset.tab.gz'},
            'geneset_json': json_object.id,
            'source': 'UCSC'
        }

        geneset.status = Data.STATUS_DONE
        geneset.save()

        with open(os.path.join(self.data_dir, str(geneset.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Generate gene set. Gene set was created '
                         'with the generate_geneset django-admin command.')

        logger.info(__('Created Gene set object: (id={})', geneset.id))

    def handle(self, *args, **options):
        """Command handle."""
        if options['rseed']:
            random.seed(42)
        for _ in range(options['n_geneset']):
            self.create_geneset()
