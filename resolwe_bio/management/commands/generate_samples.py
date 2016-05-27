import gzip
import json
import os
import random
import shutil
import string
import zipfile
import datetime

from django.core.management.base import BaseCommand
from django.contrib.auth.models import User
from django.conf import settings
from django.utils import timezone
from resolwe.flow.models import Data, DescriptorSchema, Process, Storage
from resolwe_bio.models import Sample
from resolwe_bio.tools.utils import escape_mongokey, gzopen


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
        return ''.join(random.choice(string.lowercase) for _ in range(length))

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

    def set_user(self):
        users = User.objects.filter(is_superuser=True).order_by('date_joined')
        if not users.exists():
            self.stderr.write("Admin does not exist: create a superuser")
            exit(1)
        return users.first()

    @staticmethod
    def generate_expressions(gene_ids, path):
        """Generate random expression data"""
        genes = {}
        with gzip.open(os.path.join(path, 'expressions.tab.gz'), 'wb') as f:
            f.write('{}\t{}\n'.format('Gene', 'Expression'))
            with gzopen(gene_ids) as gene_ids:
                all_genes = [line.strip() for line in gene_ids]
                for gene in all_genes:
                    expression = random.gammavariate(1, 100)
                    f.write('{}\t{}\n'.format(gene, expression))
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

    @staticmethod
    def generate_sample_desciptor(organism_name):
        annotator = ['Billie Joe Armstrong', 'Dexter Holland', 'Mark Hoppus']
        growth_protocol = ('One vial of cells was thawed out and placed in a T75 '
                           'flask in DMEM/F12 supplemented with 10% FBS and fed '
                           'every 2 days. Cells were grown to ~50% confluence '
                           '(~4-5 days), after which they were trypsinized.')
        library_prep = ('Total RNAs were extracted using the RNeasy Mini Kit (Qiagen). '
                        'RNA sequencing libraries were constructed under the standard '
                        'protocol of Illumina TruSeq RNA Prep Kit.')
        if organism_name.startswith('Dd'):
            annotation = {'geo': {'organism': 'Dictyostelium discoideum',
                                  'annotator': random.choice(annotator),
                                  'strain': random.choice(['AX2', 'AX4']),
                                  'genotype': random.choice(['wildtype', 'tirA-']),
                                  'experiment_type': 'RNA-Seq',
                                  'molecule': 'polyA RNA',
                                  'description': 'Expression profiling'},
                          'protocols': {
                            'fragmentation_method': 'sonication',
                            'growth_protocol': ('Cells were grown in nutrient media (HL-5) '
                                                'to mid-log phase prior to collection for development.'),
                            'treatment_protocol': 'Development',
                            'library_prep': ('TriZol (Life Sciences) -- Cells were scraped '
                                             '(filter) or pelleted (suspension) and disrupted '
                                             'in TriZol. Total RNA was extracted by phenol/chloroform '
                                             'as per manufacturer\'s instructions. Ribosomal RNA was '
                                             'depleted using the InDA-C method (Ovation, Nugen) reagents, '
                                             'following manufacturer\'s instructions.')}}
        if organism_name.startswith('Hs'):
            annotation = {'geo': {'organism': 'Homo sapiens',
                                  'annotator': random.choice(annotator),
                                  'source': 'brain',
                                  'cell_type': random.choice(['neuron', 'glia']),
                                  'genotype': random.choice(['wildtype', 'CDK6-/-']),
                                  'experiment_type': 'RNA-Seq',
                                  'molecule': 'polyA RNA',
                                  'optional_char': ['Cell line: DBTRG-05MG'],
                                  'description': "Human brain study. Expression profiling."},
                          'protocols': {
                            'fragmentation_method': 'sonication',
                            'growth_protocol': growth_protocol,
                            'library_prep': library_prep}}
        if organism_name.startswith('Mm'):
            annotation = {'geo': {'organism': 'Mus musculus',
                                  'annotator': random.choice(annotator),
                                  'source': 'brain',
                                  'genotype': random.choice(['wildtype', 'CDK6-/-']),
                                  'cell_type': random.choice(['neuron', 'microglia']),
                                  'experiment_type': 'RNA-Seq',
                                  'molecule': 'polyA RNA',
                                  'optional_char': ['Cell line: EOC 2'],
                                  'description': "Mouse brain study. Expression profiling."},
                          'protocols': {
                            'fragmentation_method': 'sonication',
                            'growth_protocol': growth_protocol,
                            'library_prep': library_prep}}
        return annotation

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
            name=self.set_name(),
            started=started,
            finished=started + datetime.timedelta(minutes=45),
            descriptor_schema=DescriptorSchema.objects.get(slug='reads'),
            descriptor=self.generate_reads_descriptor(),
            status='OK',
            process=Process.objects.get(slug='upload-fastq-single'),
            contributor=self.set_user(),
            input={'src': {'file': os.path.basename(reads)}},
            output={
                'fastq': {'file': os.path.basename(reads)},
                'fastqc_url': {
                    'url': 'fastqc/{}_fastqc/fastqc_report.html'.format(reads_name),
                    'refs': ['fastqc/{}_fastqc'.format(reads_name)], 'name': 'View'},
                'fastqc_archive': {'file': '{}_fastqc.zip'.format(reads_name)},
                'number': 13,
                'bases': "101"})

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

        # Create stdout file
        with open(os.path.join(data_dir, str(d.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload NGS reads. Sample was created with the generate_samples django-admin command.')

        # Get sample collection
        sample = Sample.objects.filter(data=d)[0]

        # Upload bam file
        bam = Data.objects.create(
            name='Mapping',
            process=Process.objects.get(slug='upload-bam-indexed'),
            contributor=self.set_user(),
            status='OK',
            input={
                'src': {'file': 'alignment_position_sorted.bam'},
                'src2': {'file': 'alignment_position_sorted.bam.bai'}},
            output={
                'bam': {'file': 'alignment_position_sorted.bam'},
                'bai': {'file': 'alignment_position_sorted.bam.bai'}})

        os.mkdir(os.path.join(data_dir, str(bam.id)))
        shutil.copy(bam_mapping, os.path.join(data_dir, str(bam.id)))
        shutil.copy(bai, os.path.join(data_dir, str(bam.id)))

        with open(os.path.join(data_dir, str(bam.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload BAM and BAM index (BAI) files. Sample '
                         'was created with the generate_samples django-admin command.')

        Sample.objects.filter(data=bam).delete()
        sample.data.add(bam)

        # Create expressios
        exp = Data.objects.create(
            name='Expression',
            process=Process.objects.get(slug='upload-expression'),
            contributor=self.set_user(),
            status='OK',
            input={'exp': {'file': 'expressions.tab.gz'}, 'exp_type': 'FPKM'})

        os.mkdir(os.path.join(data_dir, str(exp.id)))

        if d.name.startswith('Dd'):
            self.generate_expressions(dicty_genes, os.path.join(data_dir, str(exp.id)))
        if d.name.startswith('Hs'):
            self.generate_expressions(human_genes, os.path.join(data_dir, str(exp.id)))
        if d.name.startswith('Mm'):
            self.generate_expressions(mouse_genes, os.path.join(data_dir, str(exp.id)))

        json_object = Storage.objects.create(
            json=json.load(open(os.path.join(data_dir, str(exp.id), 'expressions.json'))),
            contributor=self.set_user(),
            data=exp)

        exp.output = {
                'exp': {'file': 'expressions.tab.gz'},
                'exp_type': 'FPKM',
                'exp_json': json_object.id
        }

        exp.save()

        Sample.objects.filter(data=exp).delete()
        sample.data.add(exp)

        with open(os.path.join(data_dir, str(exp.id), 'stdout.txt'), 'w') as stdout:
            stdout.write('Upload gene expressions. Sample was created '
                         'with the generate_samples django-admin command.')

        # Annotate Sample Collection
        if annotated:
            sample.descriptor = self.generate_sample_desciptor(d.name)
            sample.presample = False
            sample.save()
            print "Created sample: {} (id={})".format(sample.name, sample.id)
        else:
            print "Created presample: {} (id={})".format(sample.name, sample.id)

    def handle(self, *args, **options):
        if options['rseed']:
            random.seed(42)
        for i in range(options['n_samples']):
            self.create_data(annotated=True)
        for i in range(options['n_presamples']):
            self.create_data(annotated=False)
