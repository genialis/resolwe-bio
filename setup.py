import os.path

import setuptools

# Get long description from README.
with open("README.rst", "r") as fh:
    long_description = fh.read()

# Get package metadata from '__about__.py' file.
about = {}
base_dir = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(base_dir, "resolwe_bio", "__about__.py"), "r") as fh:
    exec(fh.read(), about)

setuptools.setup(
    name=about["__title__"],
    use_scm_version=True,
    description=about["__summary__"],
    long_description=long_description,
    long_description_content_type="text/x-rst",
    author=about["__author__"],
    author_email=about["__email__"],
    url=about["__url__"],
    license=about["__license__"],
    # Exclude tests from built/installed package.
    packages=setuptools.find_packages(
        exclude=["tests", "tests.*", "*.tests", "*.tests.*"]
    ),
    package_data={
        "resolwe_bio": [
            "descriptors/*.yml",
            "fixtures/*.yaml",
            "kb/migrations/*.sql",
            "migrations/*.sql",
            "processes/**/*.yml",
            "processes/**/*.py",
            "tools/*.py",
            "tools/*.R",
            "tools/*.sh",
        ]
    },
    python_requires=">=3.10, <3.12",
    install_requires=(
        "Django~=4.2",
        "djangorestframework~=3.14.0",
        "django-filter~=23.1",
        # XXX: Required due to issue https://github.com/pypa/pip/issues/4905.
        "resolwe >=38.0a1, ==38.*",
        "wrapt~=1.15.0",
    ),
    extras_require={
        "docs": [
            "daphne",
            "Sphinx~=6.1.3",
            "sphinx-rtd-theme==1.3.0",
            "pyasn1>=0.4.8",
        ],
        "package": ["twine", "wheel"],
        "test": [
            "black==23.12.0",
            "daphne",
            "flake8>=6.0.0",
            "isort>=5.12.0",
            "colorama",
            "pydocstyle~=6.3.0",
            "tblib>=1.7.0",
            "check-manifest",
            "setuptools_scm",
            "twine",
            "six==1.16",
            "build==0.10.0",
        ],
    },
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Web Environment",
        "Framework :: Django",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    keywords="bioinformatics resolwe bio pipelines dataflow django",
)
