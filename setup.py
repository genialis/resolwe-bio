#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from setuptools import setup


NAME = 'Resolwe Bioinformatics'
VERSION = '0.3.0'
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
        dependency_links=(
            "git+https://github.com/genialis/resolwe.git@9d6911e5834b5d959f80ec20906813903a35c084#egg=resolwe-0.9",
        ),
        install_requires=(
            "resolwe>=0.9",
        ),
        extras_require = {
            'docs':  ['sphinx>=1.3.2'],
            'package': [
                'twine',
                'wheel',
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
