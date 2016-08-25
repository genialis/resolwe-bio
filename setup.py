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
            'processes/**/*.yml',
            'tools/*.py',
            'tools/*.R',
        ]
    },
    zip_safe=False,
    dependency_links=(
        "git+https://github.com/genialis/resolwe.git@master#egg=resolwe-1.3.99",
    ),
    install_requires=(
        'resolwe>=1.3.1',
    ),
    extras_require = {
        'docs':  [
            'sphinx>=1.3.2',
            'sphinx_rtd_theme',
        ],
        'package': [
            'twine',
            'wheel',
        ],
        'test': [
            'check-manifest',
            'coverage>=4.2',
            'pycodestyle>=2.0.0',
            'pydocstyle>=1.0.0',
            'pylint>=1.6.4',
            'readme_renderer',
        ],
    },

    test_suite='resolwe_bio.tests',

    classifiers=[
        'Development Status :: 4 - Beta',

        'Environment :: Web Environment',
        'Framework :: Django',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Software Development :: Libraries :: Python Modules',

        'License :: OSI Approved :: Apache Software License',

        'Operating System :: OS Independent',

        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    keywords='bioinformatics resolwe bio pipelines dataflow django',
)
