from __future__ import absolute_import, division, print_function, unicode_literals

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
