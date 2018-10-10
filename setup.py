#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Use codecs' open for a consistent encoding
from codecs import open
from os import path
from setuptools import find_packages, setup


base_dir = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(base_dir, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

# Get package metadata from 'resolwe.__about__.py' file
about = {}
with open(path.join(base_dir, 'resolwe_bio', '__about__.py'), encoding='utf-8') as f:
    exec(f.read(), about)

setup(
    name=about['__title__'],

    version=about['__version__'],

    description=about['__summary__'],
    long_description=long_description,

    url=about['__url__'],

    author=about['__author__'],
    author_email=about['__email__'],

    license=about['__license__'],

    # exclude tests from built/installed package
    packages=find_packages(exclude=['tests', 'tests.*', '*.tests', '*.tests.*']),
    package_data={
        'resolwe_bio': [
            'descriptors/*.yml',
            'fixtures/*.yaml',
            'processes/**/*.yml',
            'tools/*.py',
            'tools/*.R',
            'tools/*.sh',
        ]
    },
    zip_safe=False,
    install_requires=(
        'Django~=1.11.0',
        # XXX: Remove django-autoslug after all migrations that import
        # it are deleted
        'django-autoslug==1.9.3',
        'djangorestframework~=3.7.0',
        'djangorestframework-filters~=0.10.0',
        'elasticsearch-dsl~=5.4.0',
        # XXX: Required due to issue https://github.com/pypa/pip/issues/4905.
        'resolwe >=13.0a1, ==13.*',
        'wrapt>=1.10.8',
        # XXX: djangorestframework-filters has too open requirement for
        # django-filter and doesn't work with the latest version, so we
        # have to pin it
        'django-filter~=1.0.0',
    ),
    python_requires='>=3.6, <3.7',
    extras_require = {
        'docs':  [
            # XXX: Temporarily pin Sphinx to version 1.5.x since 1.6 doesn't work with our custom
            # page template
            'Sphinx~=1.5.6',
            'sphinx_rtd_theme',
        ],
        'package': [
            'twine',
            'wheel',
        ],
        'test': [
            'check-manifest',
            # pycodestyle 2.3.0 raises false-positive for variables
            # starting with 'def'
            # https://github.com/PyCQA/pycodestyle/issues/617
            'pycodestyle~=2.2.0',
            'pydocstyle>=1.0.0',
            'pylint~=1.8.0',
            'readme_renderer',
            'tblib>=1.3.0',
        ],
    },

    test_suite='resolwe_bio.tests',

    classifiers=[
        'Development Status :: 5 - Production/Stable',

        'Environment :: Web Environment',
        'Framework :: Django',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',

        'License :: OSI Approved :: Apache Software License',

        'Operating System :: OS Independent',

        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    keywords='bioinformatics resolwe bio pipelines dataflow django',
)
