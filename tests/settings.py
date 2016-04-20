"""
Django settings for running tests for Resolwe package.

"""
from __future__ import absolute_import, division, print_function, unicode_literals

import os

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))

SECRET_KEY = 'secret'

DEBUG = True

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
)


INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.staticfiles',

    'rest_framework',
    'guardian',
    'mathfilters',
    'versionfield',

    'resolwe',
    'resolwe.permissions',
    'resolwe.flow',
    'resolwe_bio',
)

ROOT_URLCONF = 'tests.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
            ],
        },
    },
]

AUTHENTICATION_BACKENDS = (
    'django.contrib.auth.backends.ModelBackend',
    'guardian.backends.ObjectPermissionBackend',
)

ANONYMOUS_USER_ID = -1

# This is needed for runing concurrent tests on Jenkins
toxenv = os.environ.get('TOXENV', '')

# Check if PostgreSQL settings are set via environment variables
pgname = os.environ.get('RESOLWE_POSTGRESQL_NAME', 'resolwe-bio')
pguser = os.environ.get('RESOLWE_POSTGRESQL_USER', 'postgres')
pghost = os.environ.get('RESOLWE_POSTGRESQL_HOST', 'localhost')
pgport = int(os.environ.get('RESOLWE_POSTGRESQL_PORT', 5432))

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.postgresql_psycopg2',
        'NAME': pgname,
        'USER': pguser,
        'HOST': pghost,
        'PORT': pgport,
        'TEST': {
            'NAME': 'resolwe-bio_test' + toxenv
        }
    }
}

STATIC_URL = '/static/'

FLOW_EXECUTOR = {
    'NAME': 'resolwe.flow.executors.docker',
    'CONTAINER_IMAGE': 'resolwe/bio-linux8-resolwe',
    'DATA_PATH': os.path.join(PROJECT_ROOT, '.data'),
    'UPLOAD_PATH': os.path.join(PROJECT_ROOT, '.upload'),
}
# Set custom executor command if set via environment variable
if 'RESOLWE_EXECUTOR_COMMAND' in os.environ:
    FLOW_EXECUTOR['COMMAND'] = os.environ['RESOLWE_EXECUTOR_COMMAND']
FLOW_API = {
    'PERMISSIONS': 'resolwe.permissions.genesis',
}
FLOW_EXPRESSION_ENGINES = [
    'resolwe.flow.exprengines.dtlbash'
]

FLOW_DOCKER_MAPPINGS = [
    {'src': os.path.join(FLOW_EXECUTOR['DATA_PATH'], '{data_id}'),
     'dest': '/home/biolinux/data',
     'mode': 'rw'},
    {'src': FLOW_EXECUTOR['DATA_PATH'],
     'dest': '/home/biolinux/data_all',
     'mode': 'ro'},
    {'src': FLOW_EXECUTOR['UPLOAD_PATH'],
     'dest': '/home/biolinux/upload',
     'mode': 'rw'},
]

REST_FRAMEWORK = {
    'DEFAULT_AUTHENTICATION_CLASSES': (
        'rest_framework.authentication.SessionAuthentication',
    ),
    'DEFAULT_FILTER_BACKENDS': (
        'rest_framework.filters.DjangoObjectPermissionsFilter',
        'rest_framework.filters.DjangoFilterBackend',
    ),
}

FLOW_PROCESSES_FINDERS = (
    'resolwe.flow.finders.FileSystemProcessesFinder',
    'resolwe.flow.finders.AppDirectoriesFinder',
)

FLOW_PROCESSES_DIRS = (os.path.join(PROJECT_ROOT, '../resolwe_bio/tests/'),)
