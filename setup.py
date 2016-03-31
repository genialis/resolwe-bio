#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from setuptools import setup

NAME = 'resolwe-bio'
TITLE = 'Resolwe Bioinformatics'
VERSION = '1.0.0'
DESCRIPTION = "Bioinformatics pipelines for the Resolwe platform."
LONG_DESCRIPTION = open(os.path.join(os.path.dirname(__file__), 'README.rst')).read()
AUTHOR = 'Genialis'
AUTHOR_EMAIL = 'dev-team@genialis.com'
URL = 'https://github.com/genialis/resolwe-bio/'
LICENSE = 'Apache License (2.0)'

if __name__ == '__main__':
    setup(
        name=NAME,

        version=VERSION,

        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,

        url=URL,

        author=AUTHOR,
        author_email=AUTHOR_EMAIL,

        license=LICENSE,

        packages=['resolwe_bio'],
        include_package_data=True,
        zip_safe=False,
        install_requires=(
            "resolwe>=1.0.0",
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
                'django-jenkins>=0.17.0',
                'coverage>=3.7.1',
                'pep8>=1.6.2',
                'pylint>=1.4.3',
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
