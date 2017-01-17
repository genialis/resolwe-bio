""".. Ignore pydocstyle D400.

=================================
Generate Differential Expressions
=================================

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
from resolwe_bio.models import Sample

from .utils import get_descriptorschema, get_process, get_superuser, generate_sample_descriptor


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


class Command(BaseCommand):
    """Generate test differential expression data."""

    help = "Generate differential expressions"

    def __init__(self, *args, **kwargs):
        """Set command defaults."""
        super(Command, self).__init__(*args, **kwargs)

        self.data_dir = settings.FLOW_EXECUTOR['DATA_DIR']
        self.test_files_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', '..', 'tests', 'files'))

    def add_arguments(self, parser):
        """Define command arguments."""
        parser.add_argument('-n', '--n-diffexps', type=int, default=4,
                            help="Number of differential expressions to generate (default: %(default)s)")
        parser.add_argument('-g', '--group-size', type=int, default=5,
                            help="Number of samples in case/control group per DE (default: %(default)s)")
        parser.add_argument('--rseed', action='store_true', help="Use fixed random seed")

    @staticmethod
    def get_random_word(length):
        """Generate a random word."""
        return ''.join(random.choice(string.ascii_lowercase) for _ in range(length))

    def get_name(self, de_name):
        """Generate a random name."""
        return 'DE_{}_{}'.format(self.get_random_word(4), de_name)

    @staticmethod
    def generate_expressions(gene_ids, path, exp_name):
        """Generate random expression data."""
        genes = {}
        with gzip.open(os.path.join(path, exp_name), mode='wt') as f:
            # NOTE: Default line terminator is '\r\n'
            # NOTE: Python2's csv module doesn't accept a unicode string for delimeter
            csvwriter = csv.writer(f, delimiter=str('\t'), lineterminator='\n')
            csvwriter.writerow(('Gene', 'Expression'))
            with gzip.open(gene_ids, mode='rt') as gene_ids:
                all_genes = [line.strip() for line in gene_ids]
                for gene in all_genes:
                    expression = random.gammavariate(1, 100)
                    csvwriter.writerow((gene, expression))
                    genes[gene] = expression

        with open(os.path.join(path, 'expressions.json'), 'w') as json_file:
            json.dump({'genes': genes}, json_file, indent=4, sort_keys=True)

    def create_expressions(self, num, gene_ids):
        """Generate expressions."""
        expressions = []
        sample_name = 'Deseq2_{}'.format(self.get_random_word(4))

        for i in range(num):
            expression_file = 'exp_{}.tab.gz'.format(random.choice([1, 2, 3]))

            # Create expressios
            started = timezone.now()
            exp = Data.objects.create(
                name='Smpl_Ex_{}_rep{}'.format(sample_name, i + 1),
                process=get_process(slug='upload-expression'),
                contributor=get_superuser(),
                started=started,
                finished=started + datetime.timedelta(minutes=60),
                status=Data.STATUS_PROCESSING,
                input={'exp': {'file': expression_file},
                       'exp_type': 'FPKM',
                       'exp_name': 'Expression',
                       'source': 'UCSC'})

            os.mkdir(os.path.join(self.data_dir, str(exp.id)))
            self.generate_expressions(gene_ids, os.path.join(self.data_dir, str(exp.id)),
                                      expression_file)

            json_object = Storage.objects.create(
                json=json.load(open(os.path.join(self.data_dir, str(exp.id), 'expressions.json'))),
                contributor=get_superuser(),
                name='{}_storage'.format(exp.name),
                data=exp)

            exp.output = {
                'exp': {'file': expression_file},
                'exp_type': 'FPKM',
                'exp_json': json_object.id,
                'source': 'UCSC'
            }
            exp.status = Data.STATUS_DONE
            exp.save()

            sample = Sample.objects.filter(data=exp)[0]
            sample.descriptor = generate_sample_descriptor('Hs_')
            sample.annotated = True
            sample.save()

            with open(os.path.join(self.data_dir, str(exp.id), 'stdout.txt'), 'w') as stdout:
                stdout.write('Upload gene expressions. Sample was created '
                             'with the generate_diffexpr_deseq django-admin command.')

            logger.info(__('Created sample: {} (id={})', sample.name, sample.id))
            logger.info(__('\tData object: (id={})', exp.id))
            expressions.append(exp)

        return expressions

    @staticmethod
    def generate_raw_data(gene_ids, path):
        """Generate random DE data."""
        de_data = {}
        header = ['', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']

        with gzip.open(gene_ids, mode='rt') as gene_ids:
            all_genes = [line.strip() for line in gene_ids]
            de_data[''] = all_genes
            de_data['baseMean'] = [random.uniform(5, 500) for _ in all_genes]
            de_data['log2FoldChange'] = [random.uniform(-10, 10) for _ in all_genes]
            de_data['lfcSE'] = [random.uniform(0, 1) for _ in all_genes]
            de_data['stat'] = [random.uniform(-10, 10) for _ in all_genes]
            de_data['pvalue'] = [random.uniform(0, 1) for _ in all_genes]
            de_data['padj'] = [random.uniform(0, 1) for _ in all_genes]

            rows = zip(de_data[''], de_data['baseMean'], de_data['log2FoldChange'],
                       de_data['lfcSE'], de_data['stat'], de_data['pvalue'], de_data['padj'])

            with gzip.open(os.path.join(path, 'de_raw.tab.gz'), 'wt') as raw_df:
                writer = csv.writer(raw_df, delimiter=str('\t'), lineterminator='\n')
                writer.writerow(header)
                for row in rows:
                    writer.writerow(row)

        with open(os.path.join(path, 'de_json.json'), 'w') as json_file:
            de_data_std = {
                'stat': de_data['stat'],
                'logfc': de_data['log2FoldChange'],
                'pvalue': de_data['pvalue'],
                'fdr': de_data['padj'],
                'gene_id': de_data['']
            }
            json.dump(de_data_std, json_file, indent=4, sort_keys=True)

            rows = zip(de_data_std['gene_id'], de_data_std['logfc'], de_data_std['fdr'],
                       de_data_std['pvalue'], de_data_std['stat'])

            with gzip.open(os.path.join(path, 'de_file.tab.gz'), 'wt') as de_file:
                writer = csv.writer(de_file, delimiter=str('\t'), lineterminator='\n')
                writer.writerow(['gene_id', 'logfc', 'fdr', 'pvalue', 'stat'])
                for row in rows:
                    writer.writerow(row)

    def generate_diffexp_data(self, group_size):
        """Generate differential expression data."""
        de_name = 'deseq2'

        human_genes = os.path.join(self.test_files_path, 'human_genes.tab.gz')

        logger.info('---- case samples ----')
        case_samples = self.create_expressions(group_size, human_genes)

        logger.info('---- control samples ----')
        control_samples = self.create_expressions(group_size, human_genes)

        logger.info('---- upload annotation ----')

        case_input = [sample.id for sample in case_samples]
        control_input = [sample.id for sample in control_samples]

        de_descriptor = {
            'thresholds': {
                'prob_field': 'fdr',
                'logfc': 2,
                'prob': 0.05},
            'case_label': 'Case group',
            'control_label': 'Control group'
        }

        de_inputs = {
            'case': case_input,
            'control': control_input,
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
            process=get_process(slug='differentialexpression-deseq2'),
            contributor=get_superuser(),
            input=de_inputs)

        # Create data directory
        os.mkdir(os.path.join(self.data_dir, str(de_obj.id)))
        self.generate_raw_data(human_genes, os.path.join(self.data_dir, str(de_obj.id)))

        json_object = Storage.objects.create(
            json=json.load(open(os.path.join(self.data_dir, str(de_obj.id), 'de_json.json'))),
            contributor=get_superuser(),
            name='{}_storage'.format(de_obj.name),
            data=de_obj)

        os.remove(os.path.join(self.data_dir, str(de_obj.id), 'de_json.json'))

        de_obj.output = {
            'raw': {'file': 'de_raw.tab.gz'},
            'de_json': json_object.id,
            'de_file': {'file': 'de_file.tab.gz'},
            'source': 'UCSC'
        }

        de_obj.status = Data.STATUS_DONE
        de_obj.save()

        logger.info('---- new differential expression ----')
        logger.info(__('DE created with id: {}', de_obj.id))

        # Create stdout file
        with open(os.path.join(self.data_dir, str(de_obj.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Differential expression was '
                         'created with the generate_diffexpr_deseq django-admin command.')

    def handle(self, *args, **options):
        """Command handle."""
        if options['rseed']:
            random.seed(42)
        for _ in range(options['n_diffexps']):
            self.generate_diffexp_data(options['group_size'])
