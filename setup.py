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
    python_requires='>=3.6, <3.10',
    install_requires=(
        'Django~=3.1.7',
        'djangorestframework~=3.12.2',
        'django-filter~=2.4.0',
        # XXX: Required due to issue https://github.com/pypa/pip/issues/4905.
        'resolwe >=28.0a1, ==28.*',
        'wrapt~=1.12.1',
    ),
    extras_require={
        'docs': [
            'Sphinx~=3.5.3',
            'sphinx_rtd_theme',
            'pyasn1>=0.4.8',
        ],
        'package': ['twine', 'wheel'],
        'test': [
            'black',
            'flake8>=3.8.4',
            'isort>=5.7.0',
            'pydocstyle~=5.1.1',
            'tblib>=1.7.0',
            'check-manifest',
            'setuptools_scm',
            'twine',
            'six==1.15',
            # Packaging fails on Jenkins with latest version (0.3.0). Locally
            # it works fine, though.
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
    ],
    keywords='bioinformatics resolwe bio pipelines dataflow django',
)
