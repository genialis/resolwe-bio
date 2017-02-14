""".. Ignore pydocstyle D400.

================
Generate Samples
================

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
import zipfile
import datetime

from django.core.management.base import BaseCommand
from django.conf import settings
from django.utils import timezone

from resolwe.flow.models import Data, Storage
from resolwe.utils import BraceMessage as __
from resolwe_bio.models import Sample

from .utils import (get_descriptorschema, get_process, get_superuser,
                    generate_sample_descriptor, generate_reads_descriptor)


logger = logging.getLogger(__name__)  # pylint: disable=invalid-name


class Command(BaseCommand):
    """Generate test data."""

    help = "Generate test data"

    def add_arguments(self, parser):
        """Command arguments."""
        parser.add_argument('-s', '--n-samples', type=int, default=15,
                            help="Number of samples to generate (default: %(default)s)")
        parser.add_argument('-p', '--n-presamples', type=int, default=5,
                            help="Number of presamples to generate (default: %(default)s)")
        parser.add_argument('--rseed', action='store_true', help="Use fixed random seed")

    @staticmethod
    def get_random_word(length):
        """Generate a random word."""
        return ''.join(random.choice(string.ascii_lowercase) for _ in range(length))

    def set_name(self):
        """Set sample name."""
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

    def set_source(self, species):
        """Set Gene ID source."""
        if species.startswith('Dd'):
            return 'DICTYBASE'
        if species.startswith('Mm'):
            return 'UCSC'
        if species.startswith('Hs'):
            return 'UCSC'

    @staticmethod
    def generate_expressions(gene_ids, path):
        """Generate random expression data."""
        genes = {}
        with gzip.open(os.path.join(path, 'expressions.tab.gz'), mode='wt') as f:
            # NOTE: Default line terminator is '\r\n'
            # NOTE: Python2's csv module doesn't accept a unicode string for delimeter
            csvwriter = csv.writer(f, delimiter=str('\t'), lineterminator='\n')
            csvwriter.writerow(('Gene', 'Expression'))
            with gzip.open(gene_ids, mode='rt') as gene_ids:
                all_genes = [line.strip() for line in gene_ids]
                for gene in all_genes:
                    expression = round(random.gammavariate(1, 100), 2)
                    csvwriter.writerow((gene, expression))
                    genes[gene] = expression

        with open(os.path.join(path, 'expressions.json'), 'w') as json_file:
            json.dump({'genes': genes}, json_file, indent=4, sort_keys=True)

    def create_data(self, reads_name='seq_reads', annotated=False, rseed=None):
        """Generate sample data."""
        # get test data paths
        data_dir = settings.FLOW_EXECUTOR['DATA_DIR']
        test_files_path = os.path.abspath(
            os.path.join(os.path.dirname(__file__), '..', '..', 'tests', 'files'))
        reads = os.path.join(test_files_path, reads_name + '.fastq.gz')
        fastqc = os.path.join(test_files_path, reads_name + '_fastqc.zip')
        bam_mapping = os.path.join(test_files_path, 'alignment_position_sorted.bam')
        bai = os.path.join(test_files_path, 'alignment_position_sorted.bam.bai')
        dicty_genes = os.path.join(test_files_path, 'dicty_genes.tab.gz')
        human_genes = os.path.join(test_files_path, 'human_genes.tab.gz')
        mouse_genes = os.path.join(test_files_path, 'mouse_genes.tab.gz')

        # Create reads data object
        started = timezone.now()
        data_name = self.set_name()
        d = Data.objects.create(
            slug='gs-reads',
            name=data_name,
            started=started,
            finished=started + datetime.timedelta(minutes=45),
            descriptor_schema=get_descriptorschema('reads'),
            descriptor=generate_reads_descriptor(data_name, annotated=False),
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

        with zipfile.ZipFile(fastqc) as f:
            f.extractall(os.path.join(data_dir, str(d.id), 'fastqc'))

        d.output = {
            'fastq': [{'file': os.path.basename(reads)}],
            'fastqc_url': [{
                'file': 'fastqc/{}_fastqc/fastqc_report.html'.format(reads_name),
                'refs': ['fastqc/{}_fastqc'.format(reads_name)]}],
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
                   'exp_name': 'Expression',
                   'source': self.set_source(d.name)})

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
            'exp_json': json_object.id,
            'source': self.set_source(d.name)
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
            sample.descriptor = generate_sample_descriptor(d.name)
            sample.descriptor_completed = True
            sample.save()
            d.descriptor = generate_reads_descriptor(data_name, annotated=True)
            d.save()
            logger.info(__('Created annotated sample: {} (id={})', sample.name, sample.id))
        else:
            logger.info(__('Created unannotated sample: {} (id={})', sample.name, sample.id))

    def handle(self, *args, **options):
        """Command handle."""
        if options['rseed']:
            random.seed(42)
        for _ in range(options['n_samples']):
            self.create_data(annotated=True)
        for _ in range(options['n_presamples']):
            self.create_data(annotated=False)
