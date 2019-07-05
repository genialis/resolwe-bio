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
            'processes/**/*.yml',
            'processes/**/*.py',
            'tools/*.py',
            'tools/*.R',
            'tools/*.sh',
        ]
    },
    python_requires='>=3.6, <3.8',
    install_requires=(
        'Django~=2.2.0',
        'djangorestframework~=3.9.0',
        'django-filter~=2.0.0',
        'elasticsearch-dsl~=6.3.1',
        # XXX: Required due to issue https://github.com/pypa/pip/issues/4905.
        'resolwe >=18.0a1, ==18.*',
        # XXX: Temporarily pin urllib to 1.24.x, since requests 2.21.0
        # has requirement urllib3<1.25,>=1.21.1
        'urllib3~=1.24.2',
        'wrapt~=1.11.1',
    ),
    extras_require={
        'docs': [
            # XXX: Temporarily pin Sphinx to version 1.5.x since 1.6 doesn't
            # work with our custom page template.
            'Sphinx~=1.5.6',
            'sphinx_rtd_theme',
        ],
        'package': ['twine', 'wheel'],
        'test': [
            'pycodestyle~=2.5.0',
            'pydocstyle~=3.0.0',
            'pylint~=2.3.1',
            'tblib~=1.3.0',
            'check-manifest',
            'setuptools_scm',
            'twine',
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
    ],
    keywords='bioinformatics resolwe bio pipelines dataflow django',
)
