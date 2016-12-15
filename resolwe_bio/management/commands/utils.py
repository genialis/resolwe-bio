""".. Ignore pydocstyle D400.

==============================
Utility Functions for Commands
==============================

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import datetime
import gzip
import io
import os
import random
import zipfile

from django.contrib.auth import get_user_model

from resolwe.flow.models import DescriptorSchema, Process


def get_process(slug):
    """Get process by slug."""
    process_qs = Process.objects.filter(slug=slug)
    if not process_qs.exists():
        raise KeyError("Process ({}) does not exist: try running "
                       "'./manage.py register' command.".format(slug))
    return process_qs.order_by('version').last()


def get_descriptorschema(slug):
    """Get descriptor schema by slug."""
    desschema_qs = DescriptorSchema.objects.filter(slug=slug)
    if not desschema_qs.exists():
        raise KeyError("DescriptorSchema ({}) does not exist: try running "
                       "'./manage.py register' command.".format(slug))
    return desschema_qs.order_by('version').last()


def get_superuser():
    """Get a super user."""
    user_model = get_user_model()
    users = user_model.objects.filter(is_superuser=True).order_by('date_joined')
    if not users.exists():
        raise KeyError("Admin does not exist: create a superuser")
    return users.first()


def generate_sample_descriptor(organism_name):
    """Generate sample descriptor."""
    annotator = ['Billie Joe Armstrong', 'Dexter Holland', 'Mark Hoppus']
    if organism_name.startswith('Dd'):
        descriptor = {
            'sample': {
                'organism': 'Dictyostelium discoideum',
                'annotator': random.choice(annotator),
                'strain': random.choice(['AX2', 'AX4']),
                'genotype': random.choice(['wildtype', 'tirA-']),
                'molecule': 'polyA RNA',
                'description': 'Expression profiling'
            }
        }

    if organism_name.startswith('Hs'):
        descriptor = {
            'sample': {
                'organism': 'Homo sapiens',
                'annotator': random.choice(annotator),
                'source': 'Primary tumor',
                'cell_type': random.choice(['neuron', 'glia']),
                'genotype': random.choice(['wildtype', 'CDK6-/-']),
                'molecule': 'polyA RNA',
                'optional_char': ['Cell line: DBTRG-05MG'],
                'description': "Human brain study. Expression profiling."
            }
        }

    if organism_name.startswith('Mm'):
        descriptor = {
            'sample': {
                'organism': 'Mus musculus',
                'annotator': random.choice(annotator),
                'source': 'Primary tumor',
                'genotype': random.choice(['wildtype', 'CDK6-/-']),
                'cell_type': random.choice(['neuron', 'microglia']),
                'molecule': 'polyA RNA',
                'optional_char': ['Cell line: EOC 2'],
                'description': "Mouse brain study. Expression profiling."
            }
        }

    return descriptor


def generate_reads_descriptor(organism_name, presample=True):
    """Generate read data descriptor."""
    growth_protocol = ('One vial of cells was thawed out and placed in a T75 '
                       'flask in DMEM/F12 supplemented with 10% FBS and fed '
                       'every 2 days. Cells were grown to ~50% confluence '
                       '(~4-5 days), after which they were trypsinized.')
    extract_protocol = ('Total RNAs were extracted using the RNeasy Mini Kit (Qiagen).')
    library_prep = ('RNA sequencing libraries were constructed under the standard '
                    'protocol of Illumina TruSeq RNA Prep Kit.')

    reads_info = {
        'barcode': 'ATGCATGC',
        'instrument_type': 'Illumina HiSeq X',
        'barcode_removed': True,
        'seq_date': datetime.date.today().strftime('%Y-%m-%d')
    }

    if organism_name.startswith('Dd'):
        sample_descriptor = {
            'experiment_type': 'RNA-Seq',
            'protocols': {
                'fragmentation_method': 'sonication',
                'growth_protocol': ('Cells were grown in nutrient media (HL-5) '
                                    'to mid-log phase prior to collection for development.'),
                'treatment_protocol': 'Development',
                'library_prep': library_prep,
                'extract_protocol': ('TriZol (Life Sciences) -- Cells were scraped '
                                     '(filter) or pelleted (suspension) and disrupted '
                                     'in TriZol. Total RNA was extracted by phenol/chloroform '
                                     'as per manufacturer\'s instructions. Ribosomal RNA was '
                                     'depleted using the InDA-C method (Ovation, Nugen) reagents, '
                                     'following manufacturer\'s instructions.')
            }
        }

    if organism_name.startswith('Hs'):
        sample_descriptor = {
            'experiment_type': 'RNA-Seq',
            'protocols': {
                'fragmentation_method': 'sonication',
                'growth_protocol': growth_protocol,
                'library_prep': library_prep,
                'extract_protocol': extract_protocol
            }
        }

    if organism_name.startswith('Mm'):
        sample_descriptor = {
            'experiment_type': 'RNA-Seq',
            'protocols': {
                'fragmentation_method': 'sonication',
                'growth_protocol': growth_protocol,
                'library_prep': library_prep,
                'extract_protocol': extract_protocol
            }
        }

    if presample:
        return {'reads_info': reads_info}
    else:
        sample_descriptor['reads_info'] = reads_info
        return sample_descriptor


def decompress(file_name):
    """Compression-agnostic iterator.

    Iterate over files on the archive and return a tuple of
    file name, line count and file descriptor.

    Supported file formats are .tab, .gz and .zip.

    """
    if not os.path.isfile(file_name):
        raise ValueError("Can not find file '{}'".format(file_name))

    line_count = -1
    _, ext = os.path.splitext(file_name)

    _open = None
    if ext == '.tab':
        _open = open
    elif ext == '.gz':
        _open = gzip.open

    if _open:
        with _open(file_name, 'rt') as tsv_file:
            line_count = sum(1 for row in tsv_file)

        with _open(file_name, 'rt') as tsv_file:
            yield (os.path.basename(file_name), line_count, tsv_file)
    elif ext == '.zip':
        with zipfile.ZipFile(file_name) as archive:
            for entry in archive.infolist():
                if not entry.filename.endswith('.tab'):
                    continue

                if entry.filename.startswith('__MACOSX'):
                    continue

                with archive.open(entry) as tsv_file:
                    line_count = sum(1 for row in tsv_file)

                with archive.open(entry) as tsv_file:
                    yield (entry.filename, line_count, io.TextIOWrapper(tsv_file))
    else:
        raise ValueError("Unsupported file format")
