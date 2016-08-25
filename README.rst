======================
Resolwe Bioinformatics
======================

|build| |coverage| |docs| |pypi_version| |pypi_pyversions|

.. |build| image:: https://travis-ci.org/genialis/resolwe-bio.svg?branch=master
    :target: https://travis-ci.org/genialis/resolwe-bio
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

To chat with developers or ask for help, join us on Slack_.

.. _documentation: http://resolwe-bio.readthedocs.io/
.. _Slack: http://resolwe.slack.com/


Install
=======

Prerequisites
-------------

Make sure you have Python_ (2.7 or 3.4+) installed on your system. If you don't
have it yet, follow `these instructions
<https://docs.python.org/3/using/index.html>`__.

Resolwe requires PostgreSQL_ (9.4+). Many Linux distributions already include
the required version of PostgreSQL (e.g. Fedora 22+, Debian 8+, Ubuntu 15.04+)
and you can simply install it via distribution's package manager.
Otherwise, follow `these instructions
<https://wiki.postgresql.org/wiki/Detailed_installation_guides>`__.

Additionally, installing the ``psycopg2`` dependency from PyPI_ will require
having a C compiler (e.g. GCC_) as well as Python and PostgreSQL development
files installed on the system.

Note
^^^^

The preferred way to install the C compiler and Python and PostgreSQL
development files is to use your distribution's packages, if they exist. For
example, on a Fedora/RHEL-based system, that would mean installing ``gcc``,
``python-devel``/``python3-devel`` and ``postgresql-devel`` packages.

.. _Python: https://www.python.org/
.. _PostgreSQL: http://www.postgresql.org/
.. _PyPi: https://pypi.python.org/
.. _GCC: https://gcc.gnu.org/

From PyPI_
----------

.. code::

    pip install resolwe-bio

From source
-----------

.. code::

   pip install https://github.com/genialis/resolwe-bio/archive/<git-tree-ish>.tar.gz

where ``<git-tree-ish>`` can represent any commit SHA, branch name, tag name,
etc. in `Resolwe Bioinformatics' GitHub repository`_. For example, to install
the latest Resolwe Bioinformatics from the ``master`` branch, use:

.. code::

   pip install https://github.com/genialis/resolwe-bio/archive/master.tar.gz

.. _`Resolwe Bioinformatics' GitHub repository`: https://github.com/genialis/resolwe-bio/


Contribute
==========

We welcome new contributors. To learn more, read Contributing_ section of the
documentation.

.. _Contributing: http://resolwe-bio.readthedocs.io/en/latest/contributing.html
