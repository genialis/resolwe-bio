import os.path
import setuptools

# Get long description from README.
with open('README.rst', 'r') as fh:
    long_description = fh.read()

# Get package metadata from '__about__.py' file.
about = {}
base_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(base_dir, 'resolwe_bio', '__about__.py'), 'r') as fh:
    exec(fh.read(), about)

setuptools.setup(
    name=about['__title__'],
    use_scm_version=True,
    description=about['__summary__'],
    long_description=long_description,
    long_description_content_type='text/x-rst',
    author=about['__author__'],
    author_email=about['__email__'],
    url=about['__url__'],
    license=about['__license__'],
    # Exclude tests from built/installed package.
    packages=setuptools.find_packages(
        exclude=['tests', 'tests.*', '*.tests', '*.tests.*']
    ),
    package_data={
        'resolwe_bio': [
            'descriptors/*.yml',
            'fixtures/*.yaml',
            'kb/migrations/*.sql',
            'migrations/*.sql',
            'processes/**/*.yml',
            'processes/**/*.py',
            'tools/*.py',
            'tools/*.R',
            'tools/*.sh',
        ]
    },
    python_requires='>=3.6, <3.11',
    install_requires=(
        'Django~=3.2.12',
        'djangorestframework~=3.13.1',
        'django-filter~=21.1',
        # XXX: Required due to issue https://github.com/pypa/pip/issues/4905.
        'resolwe >=31.0a1, ==31.*',
        'wrapt~=1.13.3',
    ),
    extras_require={
        'docs': [
            'Sphinx~=4.3.2',
            'sphinx_rtd_theme',
            'pyasn1>=0.4.8',
        ],
        'package': ['twine', 'wheel'],
        'test': [
            'black',
            'flake8>=4.0.1',
            'isort>=5.10.1',
            'pydocstyle~=6.1.1',
            'tblib>=1.7.0',
            'check-manifest',
            'setuptools_scm',
            'twine',
            'six==1.16',
            'build==0.2.1',
        ],
    },
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
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    keywords='bioinformatics resolwe bio pipelines dataflow django',
)
