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
import shutil
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

    def create_expressions(self, num):
        """Generate expressions."""
        expressions = []
        sample_name = 'Cuffdiff_{}'.format(self.get_random_word(4))

        for i in range(num):
            cuffquant_file = 'cuffquant_{}.cxb'.format(random.choice([1, 2]))

            # Create expressios
            exp = Data.objects.create(
                name='Smpl_Ex_{}_rep{}'.format(sample_name, i + 1),
                process=get_process('upload-cxb'),
                contributor=get_superuser(),
                status=Data.STATUS_PROCESSING,
                input={'src': {'file': cuffquant_file}, 'source': 'UCSC'})

            os.mkdir(os.path.join(self.data_dir, str(exp.id)))
            shutil.copy(os.path.join(self.test_files_path, cuffquant_file), os.path.join(self.data_dir, str(exp.id)))

            exp.output = {
                'cxb': {'file': cuffquant_file},
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
                             'with the generate_diffexr_cuffdiff django-admin command.')

            logger.info(__('Created sample: {} (id={})', sample.name, sample.id))
            logger.info(__('\tData object: (id={})', exp.id))
            expressions.append(exp)

        return expressions

    def create_genome_annotation(self, filename):
        """Create a genome annotation."""
        ann = Data.objects.create(
            name='Annotation_{}'.format(filename.split('.')[0]),
            process=get_process('upload-gtf'),
            contributor=get_superuser(),
            status=Data.STATUS_PROCESSING,
            input={'src': {'file': filename}, 'source': 'UCSC'})

        os.mkdir(os.path.join(self.data_dir, str(ann.id)))

        with gzip.open(os.path.join(self.test_files_path, filename), 'rb') as gzfile:
            with open(os.path.join(self.data_dir, str(ann.id), filename[:-3]), 'wb') as outfile:
                shutil.copyfileobj(gzfile, outfile)

        ann.output = {
            'gtf': {'file': filename[:-3]},
            'source': 'UCSC'
        }
        ann.status = Data.STATUS_DONE
        ann.save()

        with open(os.path.join(self.data_dir, str(ann.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload genome annotation with the '
                         'generate_diffexpr_cuffdiff django-admin command.')

        logger.info(__('Genome annotation created: {} (id={})', filename, ann.id))

        return ann

    @staticmethod
    def generate_raw_data(gene_ids, path):
        """Generate random DE data."""
        de_data = {}
        header = ['test_id', 'gene_id', 'gene', 'locus', 'sample_1',
                  'sample_2', 'status', 'value_1', 'value_2',
                  'log2(fold_change)', 'test_stat', 'p_value', 'q_value',
                  'significant']

        with gzip.open(gene_ids, mode='rt') as gene_ids:
            all_genes = [line.strip() for line in gene_ids]
            n_of_genes = len(all_genes)
            de_data['test_id'] = all_genes
            de_data['gene_id'] = all_genes
            de_data['gene'] = all_genes
            de_data['locus'] = ['chr20:463337-524482'] * n_of_genes
            de_data['sample_1'] = ['control'] * n_of_genes
            de_data['sample_2'] = ['case'] * n_of_genes
            de_data['status'] = ['OK'] * n_of_genes
            de_data['value_1'] = [random.gammavariate(1, 100) for _ in all_genes]
            de_data['value_2'] = [random.gammavariate(1, 100) for _ in all_genes]
            de_data['log2(fold_change)'] = [random.uniform(-10, 10) for _ in all_genes]
            de_data['test_stat'] = [random.uniform(-3, 3) for _ in all_genes]
            de_data['p_value'] = [random.uniform(0, 1) for _ in all_genes]
            de_data['q_value'] = [random.uniform(0, 1) for _ in all_genes]
            de_data['significant'] = [random.choice(['yes', 'no']) for _ in all_genes]

            rows = zip(de_data['test_id'], de_data['gene_id'], de_data['gene'],
                       de_data['locus'], de_data['sample_1'], de_data['sample_2'],
                       de_data['status'], de_data['value_1'], de_data['value_2'],
                       de_data['log2(fold_change)'], de_data['test_stat'],
                       de_data['p_value'], de_data['q_value'], de_data['significant'])

            with gzip.open(os.path.join(path, 'de_raw.tab.gz'), 'wt') as raw_df:
                writer = csv.writer(raw_df, delimiter=str('\t'), lineterminator='\n')
                writer.writerow(header)
                for row in rows:
                    writer.writerow(row)

        with open(os.path.join(path, 'de_json.json'), 'w') as json_file:
            de_data_std = {
                'stat': de_data['test_stat'],
                'logfc': de_data['log2(fold_change)'],
                'pvalue': de_data['p_value'],
                'fdr': de_data['q_value'],
                'gene_id': de_data['gene_id']
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
        de_name = 'cuffdiff'

        human_genes = os.path.join(self.test_files_path, 'human_genes.tab.gz')

        logger.info('---- case samples ----')
        case_samples = self.create_expressions(group_size)

        logger.info('---- control samples ----')
        control_samples = self.create_expressions(group_size)

        logger.info('---- upload annotation ----')

        case_input = [sample.id for sample in case_samples]
        control_input = [sample.id for sample in control_samples]
        genome_annotation_input = self.create_genome_annotation('hg19_chr20_small.gtf.gz')

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

        # Create data directory
        os.mkdir(os.path.join(self.data_dir, str(de_obj.id)))
        self.generate_raw_data(human_genes, os.path.join(self.data_dir, str(de_obj.id)))

        with open(os.path.join(self.data_dir, str(de_obj.id), 'test.txt'), 'w') as _:
            pass

        json_object = Storage.objects.create(
            json=json.load(open(os.path.join(self.data_dir, str(de_obj.id), 'de_json.json'))),
            contributor=get_superuser(),
            name='{}_storage'.format(de_obj.name),
            data=de_obj)

        os.remove(os.path.join(self.data_dir, str(de_obj.id), 'de_json.json'))

        # TODO: reference on existing true files
        de_obj.output = {
            'raw': {'file': 'de_raw.tab.gz'},
            'de_json': json_object.id,
            'de_file': {'file': 'de_file.tab.gz'},
            'transcript_diff_exp': {'file': 'test.txt'},
            'cds_diff_exp': {'file': 'test.txt'},
            'tss_group_diff_exp': {'file': 'test.txt'},
            'cuffdiff_output': {'file': 'test.txt'},
            'source': 'UCSC'
        }

        de_obj.status = Data.STATUS_DONE
        de_obj.save()

        logger.info('---- new differential expression ----')
        logger.info(__('DE created with id: {}', de_obj.id))

        # Create stdout file
        with open(os.path.join(self.data_dir, str(de_obj.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Differential expression was '
                         'created with the generate_diffexpr_cuffdiff django-admin command.')

    def handle(self, *args, **options):
        """Command handle."""
        if options['rseed']:
            random.seed(42)
        for _ in range(options['n_diffexps']):
            self.generate_diffexp_data(options['group_size'])
