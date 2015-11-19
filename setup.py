#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from setuptools import setup


NAME = 'Resolwe Bioinformatics'
VERSION = '0.0.1'
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
        author=AUTHOR,
        author_email=AUTHOR_EMAIL,
        url=URL,
        license=LICENSE,
        packages=['resolwe_bio'],
        include_package_data=True,
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Web Environment',
            'Framework :: Django',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Apache Software License',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Software Development :: Libraries :: Python Modules',
        ],
        zip_safe=False,
        install_requires=(
            # "resolwe>=0.0.1",
        ),
        test_suite='resolwe_bio.tests',
    )
