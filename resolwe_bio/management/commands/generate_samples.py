from __future__ import absolute_import, division, print_function, unicode_literals

import gzip
import json
import logging
import os
import random
import shutil
import string
import zipfile
import datetime

from django.core.management.base import BaseCommand
from django.conf import settings
from django.utils import timezone

from resolwe.flow.models import Data, Storage
from resolwe_bio.models import Sample
from resolwe_bio.tools.utils import escape_mongokey, gzopen
from .utils import get_descriptorschema, get_process, get_superuser, generate_sample_desciptor


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


class Command(BaseCommand):

    """Generate test data"""

    help = "Generate test data"

    def add_arguments(self, parser):
        parser.add_argument('-s', '--n-samples', type=int, default=15,
                            help="Number of samples to generate (default: %(default)s)")
        parser.add_argument('-p', '--n-presamples', type=int, default=5,
                            help="Number of presamples to generate (default: %(default)s)")
        parser.add_argument('--rseed', action='store_true', help="Use fixed random seed")

    @staticmethod
    def get_random_word(length):
        return ''.join(random.choice(string.ascii_lowercase) for _ in range(length))

    def set_name(self):
        organism = random.choice(['Dictyostelium discoideum', 'Mus musculus', 'Homo sapiens'])
        replicate = random.choice(['rep1', 'rep2', 'rep3', 'rep4', 'rep5'])
        hour = random.choice(range(36))
        kit = random.choice(['RiboZero', 'Nugen'])
        group = random.choice(['treatment', 'control'])
        if organism == 'Dictyostelium discoideum':
            return 'Dd_{}_{}_hr_{}'.format(kit, replicate, hour)
        if organism == 'Mus musculus':
            return 'Mm_{}_{}_{}'.format(self.get_random_word(3), group, replicate)
        if organism == 'Homo sapiens':
            return 'Hs_{}_{}_{}'.format(self.get_random_word(3), group, replicate)

    @staticmethod
    def generate_expressions(gene_ids, path):
        """Generate random expression data"""
        genes = {}
        with gzip.open(os.path.join(path, 'expressions.tab.gz'), 'wb') as f:
            f.write('{}\t{}\n'.format('Gene', 'Expression').encode('utf-8'))
            with gzopen(gene_ids) as gene_ids:
                all_genes = [line.strip() for line in gene_ids]
                for gene in all_genes:
                    expression = random.gammavariate(1, 100)
                    f.write('{}\t{}\n'.format(gene, expression).encode('utf-8'))
                    genes[escape_mongokey(gene)] = expression

        with open(os.path.join(path, 'expressions.json'), 'w') as json_file:
            js_out = '{{"exp_json":{}}}'.format(json.dumps({'genes': genes}, separators=(',', ':')))
            json_file.write(js_out)

    @staticmethod
    def generate_reads_descriptor():
        barcodes = ['ATGCATGC', 'TACGTACG']
        instruments = ['Illumina HiSeq X', 'Illumina HiSeq 3000', 'Illumina HiSeq 3000']
        descriptor = {'barcode': random.choice(barcodes),
                      'instrument_type': random.choice(instruments),
                      'barcode_removed': True}

        seq_date = datetime.date.today().strftime('%Y-%m-%d')
        descriptor['seq_date'] = seq_date

        return descriptor

    def create_data(self, reads_name='seq_reads', annotated=False, rseed=None):
        # get test data paths
        data_dir = settings.FLOW_EXECUTOR['DATA_DIR']
        test_files_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', '..', 'tests', 'processes', 'files'))
        reads = os.path.join(test_files_path, reads_name + '.fastq.gz')
        fastqc = os.path.join(test_files_path, reads_name + '_fastqc.zip')
        bam_mapping = os.path.join(test_files_path, 'alignment_position_sorted.bam')
        bai = os.path.join(test_files_path, 'alignment_position_sorted.bam.bai')
        dicty_genes = os.path.join(test_files_path, 'dicty_genes.tab.gz')
        human_genes = os.path.join(test_files_path, 'human_genes.tab.gz')
        mouse_genes = os.path.join(test_files_path, 'mouse_genes.tab.gz')

        # Create reads data object
        started = timezone.now()
        d = Data.objects.create(
            slug='gs-reads',
            name=self.set_name(),
            started=started,
            finished=started + datetime.timedelta(minutes=45),
            descriptor_schema=get_descriptorschema('reads'),
            descriptor=self.generate_reads_descriptor(),
            status=Data.STATUS_PROCESSING,
            process=get_process('upload-fastq-single'),
            contributor=get_superuser(),
            input={'src': [{'file': os.path.basename(reads)}]})

        # Create data directory and copy reads files into it
        os.mkdir(os.path.join(data_dir, str(d.id)))
        shutil.copy(reads, os.path.join(data_dir, str(d.id)))

        # Attach FastQC data to reads file
        os.mkdir(os.path.join(data_dir, str(d.id), 'fastqc'))
        shutil.copy(fastqc, os.path.join(data_dir, str(d.id)))

        with zipfile.ZipFile(fastqc) as zf:
            zf.extractall(os.path.join(data_dir, str(d.id), 'fastqc'))

        old_fastqc_path = os.path.join(data_dir, str(d.id), 'fastqc', reads_name + '_fastqc', 'fastqc_report.html')
        new_fastqc_path = os.path.join(data_dir, str(d.id), 'fastqc', reads_name + '_fastqc.html')
        shutil.copy(old_fastqc_path, new_fastqc_path)

        d.output = {
            'fastq': [{'file': os.path.basename(reads)}],
            'fastqc_url': [{
                'url': 'fastqc/{}_fastqc/fastqc_report.html'.format(reads_name),
                'refs': ['fastqc/{}_fastqc'.format(reads_name)], 'name': 'View'}],
            'fastqc_archive': [{'file': '{}_fastqc.zip'.format(reads_name)}]}
        d.status = Data.STATUS_DONE
        d.save()

        # Create stdout file
        with open(os.path.join(data_dir, str(d.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload NGS reads. Sample was created with the generate_samples django-admin command.')

        # Get sample collection
        sample = Sample.objects.filter(data=d)[0]

        # Upload bam file
        bam = Data.objects.create(
            name='Mapping',
            started=started,
            finished=started + datetime.timedelta(minutes=50),
            process=get_process('upload-bam-indexed'),
            contributor=get_superuser(),
            status=Data.STATUS_PROCESSING,
            input={
                'src': {'file': 'alignment_position_sorted.bam'},
                'src2': {'file': 'alignment_position_sorted.bam.bai'}})

        os.mkdir(os.path.join(data_dir, str(bam.id)))
        shutil.copy(bam_mapping, os.path.join(data_dir, str(bam.id)))
        shutil.copy(bai, os.path.join(data_dir, str(bam.id)))

        bam.output = {
            'bam': {'file': 'alignment_position_sorted.bam'},
            'bai': {'file': 'alignment_position_sorted.bam.bai'}}
        bam.status = Data.STATUS_DONE
        bam.save()

        with open(os.path.join(data_dir, str(bam.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload BAM and BAM index (BAI) files. Sample '
                         'was created with the generate_samples django-admin command.')

        Sample.objects.filter(data=bam).delete()
        sample.data.add(bam)

        # Create expressios
        exp = Data.objects.create(
            name='Expression',
            process=get_process(slug='upload-expression'),
            contributor=get_superuser(),
            started=started,
            finished=started + datetime.timedelta(minutes=60),
            status=Data.STATUS_PROCESSING,
            input={'exp': {'file': 'expressions.tab.gz'},
                   'exp_type': 'FPKM',
                   'exp_name': 'Expression'})

        os.mkdir(os.path.join(data_dir, str(exp.id)))

        if d.name.startswith('Dd'):
            self.generate_expressions(dicty_genes, os.path.join(data_dir, str(exp.id)))
        if d.name.startswith('Hs'):
            self.generate_expressions(human_genes, os.path.join(data_dir, str(exp.id)))
        if d.name.startswith('Mm'):
            self.generate_expressions(mouse_genes, os.path.join(data_dir, str(exp.id)))

        json_object = Storage.objects.create(
            json=json.load(open(os.path.join(data_dir, str(exp.id), 'expressions.json'))),
            contributor=get_superuser(),
            name='{}_storage'.format(exp.name),
            data=exp)

        exp.output = {
            'exp': {'file': 'expressions.tab.gz'},
            'exp_type': 'FPKM',
            'exp_json': json_object.id
        }
        exp.status = Data.STATUS_DONE
        exp.save()

        Sample.objects.filter(data=exp).delete()
        sample.data.add(exp)

        with open(os.path.join(data_dir, str(exp.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload gene expressions. Sample was created '
                         'with the generate_samples django-admin command.')

        # Annotate Sample Collection
        if annotated:
            sample.descriptor = generate_sample_desciptor(d.name)
            sample.presample = False
            sample.save()
            logger.info("Created sample: {} (id={})".format(sample.name, sample.id))
        else:
            logger.info("Created presample: {} (id={})".format(sample.name, sample.id))

    def handle(self, *args, **options):
        if options['rseed']:
            random.seed(42)
        for i in range(options['n_samples']):
            self.create_data(annotated=True)
        for i in range(options['n_presamples']):
            self.create_data(annotated=False)
