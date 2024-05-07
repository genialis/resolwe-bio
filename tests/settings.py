"""
Django settings for running tests for Resolwe package.

"""

import os
import re
import sys
from distutils.util import strtobool  # pylint: disable=import-error,no-name-in-module

PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))

SECRET_KEY = "secret"

# TODO: Remove this setting completely and only set it in the tests that require it.
RESOLWE_HOST_URL = "https://dummy.host.local"

DEBUG = True

ALLOWED_HOSTS = ["*"]

MIDDLEWARE_CLASSES = (
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
)


INSTALLED_APPS = (
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.staticfiles",
    "channels",
    "rest_framework",
    "versionfield",
    "resolwe",
    "resolwe.permissions",
    "resolwe.flow",
    "resolwe.observers",
    "resolwe.storage",
    "resolwe.toolkit",
    "resolwe.test_helpers",
    "resolwe_bio",
    "resolwe_bio.kb",
    "resolwe_bio.variants",
)

ROOT_URLCONF = "tests.urls"

TEST_RUNNER = "resolwe.test_helpers.test_runner.ResolweRunner"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
            ],
        },
    },
]

AUTHENTICATION_BACKENDS = (
    "django.contrib.auth.backends.ModelBackend",
    "resolwe.permissions.permissions.ResolwePermissionBackend",
)

DEFAULT_AUTO_FIELD = "django.db.models.AutoField"

ANONYMOUS_USER_NAME = "public"

# Check if PostgreSQL settings are set via environment variables
pgname = os.environ.get("RESOLWE_POSTGRESQL_NAME", "resolwe-bio")
pguser = os.environ.get("RESOLWE_POSTGRESQL_USER", "resolwe")
pgpass = os.environ.get("RESOLWE_POSTGRESQL_PASS", "resolwe")
pghost = os.environ.get("RESOLWE_POSTGRESQL_HOST", "localhost")
pgport = int(os.environ.get("RESOLWE_POSTGRESQL_PORT", 55433))

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.postgresql",
        "NAME": pgname,
        "USER": pguser,
        "PASSWORD": pgpass,
        "HOST": pghost,
        "PORT": pgport,
    }
}

STATIC_URL = "/static/"

REDIS_CONNECTION = {
    "host": "localhost",
    "port": int(os.environ.get("RESOLWE_REDIS_PORT", 56380)),
    "db": int(os.environ.get("RESOLWE_REDIS_DATABASE", 0)),
    "protocol": (os.environ.get("RESOLWE_REDIS_PROTOCOL", "redis")),
}
REDIS_CONNECTION_STRING = "{protocol}://{host}:{port}/{db}".format(**REDIS_CONNECTION)

LISTENER_CONNECTION = {
    # Keys in the hosts dictionary are workload connector names. Currently
    # supported are 'local', 'kubertenes', 'celery' and 'slurm'.
    "hosts": {"local": "172.17.0.1"},
    "port": int(os.environ.get("RESOLWE_LISTENER_SERVICE_PORT", 53893)),
    "min_port": 50000,
    "max_port": 60000,
    "protocol": "tcp",
}

# The IP address where listener is available from the communication container.
# The setting is a dictionary where key is the name of the workload connector.
COMMUNICATION_CONTAINER_LISTENER_CONNECTION = {"local": "172.17.0.1"}

# Settings in OSX/Windows are different since Docker runs in a virtual machine.
if sys.platform == "darwin":
    LISTENER_CONNECTION["hosts"]["local"] = "0.0.0.0"
    COMMUNICATION_CONTAINER_LISTENER_CONNECTION = {"local": "127.0.0.1"}

FLOW_EXECUTOR = {
    "NAME": "resolwe.flow.executors.docker",
    # XXX: Change to a stable resolwe image when it will include all the required tools
    "CONTAINER_IMAGE": "resolwe/bio-linux8-resolwe-preview",
    "CONTAINER_NAME_PREFIX": "resolwebio",
    "REDIS_CONNECTION": REDIS_CONNECTION,
    "LISTENER_CONNECTION": LISTENER_CONNECTION,
}
# Set custom executor command if set via environment variable
if "RESOLWE_DOCKER_COMMAND" in os.environ:
    FLOW_DOCKER_COMMAND = os.environ["RESOLWE_DOCKER_COMMAND"]
FLOW_API = {
    "PERMISSIONS": "resolwe.permissions.permissions",
}
FLOW_EXPRESSION_ENGINES = [
    {
        "ENGINE": "resolwe.flow.expression_engines.jinja",
        "CUSTOM_FILTERS": [
            "resolwe_bio.expression_filters.sample",
            "resolwe_bio.expression_filters.relation",
        ],
    },
]
FLOW_EXECUTION_ENGINES = [
    "resolwe.flow.execution_engines.bash",
    "resolwe.flow.execution_engines.workflow",
    "resolwe.flow.execution_engines.python",
]

# Check if any Manager settings are set via environment variables
manager_prefix = os.environ.get("RESOLWE_MANAGER_REDIS_PREFIX", "resolwe-bio.manager")
# Ensure Manager channel prefix is a valid Django Channels name.
manager_prefix = re.sub("[^0-9a-zA-Z.-]", "-", manager_prefix)
FLOW_MANAGER = {
    "NAME": "resolwe.flow.managers.workload_connectors.local",
    "REDIS_PREFIX": manager_prefix,
    "REDIS_CONNECTION": REDIS_CONNECTION,
}

FLOW_PROCESS_MAX_CORES = 1
FLOW_PROCESS_MAX_MEMORY = 10240

# Don't pull Docker images if set via the environment variable.
FLOW_DOCKER_DONT_PULL = strtobool(os.environ.get("RESOLWE_DOCKER_DONT_PULL", "0"))

# Ensure all container images follow a specific format.
FLOW_CONTAINER_VALIDATE_IMAGE = r".+:(?!latest)"

REST_FRAMEWORK = {
    "DEFAULT_AUTHENTICATION_CLASSES": (
        "rest_framework.authentication.SessionAuthentication",
    ),
    "DEFAULT_FILTER_BACKENDS": (
        "resolwe.permissions.filters.ResolwePermissionsFilter",
        "django_filters.rest_framework.DjangoFilterBackend",
        "resolwe.flow.filters.OrderingFilter",
    ),
    # Python<3.7 cannot parse iso-8601 formatted datetimes with tz-info form
    # "+01:00" (DRF default). It can only parse "+0100" form, so we need to
    # modify this setting. This will be fixed in Python3.7, where "+01:00" can
    # be parsed by ``datetime.datetime.strptime`` syntax.
    # For more, check "%z" syntax description in:
    # https://docs.python.org/3.7/library/datetime.html#strftime-and-strptime-behavior
    "DATETIME_FORMAT": "%Y-%m-%dT%H:%M:%S.%f%z",
}

# Time

USE_TZ = True

TIME_ZONE = "UTC"

# Django does not support parsing of 'iso-8601' formated datetimes by default.
# Since Django-filters uses Django forms for parsing, we need to modify Django
# setting ``DATETIME_INPUT_FORMATS`` to support 'iso-8601' format.
# https://docs.djangoproject.com/en/1.11/ref/settings/#datetime-input-formats
DATETIME_INPUT_FORMATS = (
    # These are already given Django defaults:
    "%Y-%m-%d %H:%M:%S",  # '2006-10-25 14:30:59'
    "%Y-%m-%d %H:%M:%S.%f",  # '2006-10-25 14:30:59.000200'
    "%Y-%m-%d %H:%M",  # '2006-10-25 14:30'
    "%Y-%m-%d",  # '2006-10-25'
    # These are iso-8601 formatted:
    "%Y-%m-%dT%H:%M:%S.%f%z",  # '2006-10-25T14:30:59.000200+0200' or '2006-10-25T14:30:59.000200+02:00' (Python>=3.7)
    "%Y-%m-%dT%H:%M:%S.%fZ",  # '2006-10-25T14:30:59.000200Z'
    "%Y-%m-%dT%H:%M:%S.%f",  # '2006-10-25T14:30:59.000200'
    "%Y-%m-%dT%H:%M:%SZ",  # '2006-10-25T14:30:59Z'
    "%Y-%m-%dT%H:%M:%S",  # '2006-10-25T14:30:59'
    "%Y-%m-%dT%H:%M",  # '2006-10-25T14:30'
)

FLOW_PROCESSES_FINDERS = (
    "resolwe.flow.finders.FileSystemProcessesFinder",
    "resolwe.flow.finders.AppDirectoriesFinder",
)

FLOW_PROCESSES_RUNTIMES = (
    "resolwe.process.runtime.Process",
    "resolwe_bio.process.runtime.ProcessBio",
)

FLOW_PROCESSES_DIRS = (os.path.join(PROJECT_ROOT, "../resolwe_bio/processes/"),)

# Do not skip tests that fail on Docker executor if this is set via environment
# variable
if os.environ.get("RESOLWEBIO_TESTS_SKIP_DOCKER_FAILURES", "").lower() in [
    "no",
    "false",
]:
    TESTS_SKIP_DOCKER_FAILURES = False

# Testing.

TEST_RUNNER = "resolwe.test_helpers.test_runner.ResolweRunner"
TEST_PROCESS_REQUIRE_TAGS = True
# Don't profile unless set via the environment variable.
TEST_PROCESS_PROFILE = strtobool(os.environ.get("RESOLWE_TEST_PROCESS_PROFILE", "0"))

# Channels.

ASGI_APPLICATION = "tests.routing.channel_routing"

CHANNEL_LAYERS = {
    "default": {
        "BACKEND": "channels_redis.core.RedisChannelLayer",
        "CONFIG": {
            "hosts": [REDIS_CONNECTION_STRING],
            "expiry": 3600,
        },
    },
}


# Logging.

# Set RESOLWEBIO_LOG_FILE environment variable to a file path to enable logging
# debugging messages to to a file.
log_file_path = os.environ.get(
    "RESOLWEBIO_LOG_FILE", os.devnull
)  # pylint: disable=invalid-name

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "standard": {
            "format": "%(asctime)s - %(levelname)s - %(name)s[%(process)s]: %(message)s",
        },
    },
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "level": "WARNING",
            "formatter": "standard",
        },
        "file": {
            "class": "logging.handlers.RotatingFileHandler",
            "filename": log_file_path,
            "formatter": "standard",
            "maxBytes": 1024 * 1024 * 10,  # 10 MB
        },
    },
    "loggers": {
        "": {
            "handlers": ["file"],
            "level": "DEBUG",
        },
    },
}
