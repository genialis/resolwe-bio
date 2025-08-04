======================
Resolwe Bioinformatics
======================

|build| |coverage| |docs| |pypi_version| |pypi_pyversions|

.. |build| image:: https://github.com/genialis/resolwe-bio/actions/workflows/ci.yml/badge.svg?branch=master
    :target: https://github.com/genialis/resolwe-bio/actions?query=branch%3Amaster
    :alt: Build Status

.. |coverage| image:: https://img.shields.io/codecov/c/github/genialis/resolwe-bio/master.svg
    :target: http://codecov.io/github/genialis/resolwe-bio?branch=master
    :alt: Coverage Status

.. |docs| image:: https://readthedocs.org/projects/resolwe-bio/badge/?version=latest
    :target: http://resolwe-bio.readthedocs.io/
    :alt: Documentation Status

.. |pypi_version| image:: https://img.shields.io/pypi/v/resolwe-bio.svg
    :target: https://pypi.python.org/pypi/resolwe-bio
    :alt: Version on PyPI

.. |pypi_pyversions| image:: https://img.shields.io/pypi/pyversions/resolwe-bio.svg
    :target: https://pypi.python.org/pypi/resolwe-bio
    :alt: Supported Python versions

.. |pypi_downloads| image:: https://img.shields.io/pypi/dm/resolwe-bio.svg
    :target: https://pypi.python.org/pypi/resolwe-bio
    :alt: Number of downloads from PyPI

Bioinformatics pipelines for the Resolwe_ dataflow package for `Django
framework`_.

.. _Resolwe: https://github.com/genialis/resolwe
.. _Django framework: https://www.djangoproject.com/


Docs & Help
===========

Read about getting started and how to write `processes` in the documentation_.

.. _documentation: http://resolwe-bio.readthedocs.io/


Install
=======

Prerequisites
-------------

Make sure you have Python_ 3.10 or later installed on your system. If you don't have it
yet, follow `these instructions
<https://docs.python.org/3/using/index.html>`__.

.. _Python: https://www.python.org/


Using PyPI_
-----------

.. code::

    pip install resolwe-bio

To install a pre-release, use:

.. code::

    pip install --pre resolwe-bio

.. _PyPi: https://pypi.python.org/


Using source on GitHub_
-----------------------

.. code::

   pip install --pre https://github.com/genialis/resolwe-bio/archive/<git-tree-ish>.tar.gz

where ``<git-tree-ish>`` can represent any commit SHA, branch name, tag name,
etc. in `Resolwe Bioinformatics' GitHub repository`_. For example, to install
the latest Resolwe Bioinformatics from the ``master`` branch, use:

.. code::

   pip install --pre https://github.com/genialis/resolwe-bio/archive/master.tar.gz

.. _`Resolwe Bioinformatics' GitHub repository`: https://github.com/genialis/resolwe-bio/
.. _GitHub: `Resolwe Bioinformatics' GitHub repository`_


Contribute
==========

We welcome new contributors. To learn more, read Contributing_ section of the
documentation.

.. _Contributing: http://resolwe-bio.readthedocs.io/en/latest/contributing.html
