from __future__ import absolute_import, division, print_function, unicode_literals

import random

from django.contrib.auth import get_user_model

from resolwe.flow.models import DescriptorSchema, Process


def get_process(slug):
    process_qs = Process.objects.filter(slug=slug)
    if not process_qs.exists():
        raise KeyError("Process ({}) does not exist: try running "
                       "'./manage.py register' command.".format(slug))
    return process_qs.order_by('version').last()


def get_descriptorschema(slug):
    desschema_qs = DescriptorSchema.objects.filter(slug=slug)
    if not desschema_qs.exists():
        raise KeyError("DescriptorSchema ({}) does not exist: try running "
                       "'./manage.py register' command.".format(slug))
    return desschema_qs.order_by('version').last()


def get_superuser():
    User = get_user_model()
    users = User.objects.filter(is_superuser=True).order_by('date_joined')
    if not users.exists():
        raise KeyError("Admin does not exist: create a superuser")
    return users.first()

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
