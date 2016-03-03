======================
Resolwe Bioinformatics
======================

Bioinformatics pipelines for the Resolwe_ dataflow package for `Django
framework`_.

.. _Resolwe: https://github.com/genialis/resolwe
.. _Django framework: https://www.djangoproject.com/

Docs & Help
===========

Read detailed description in the documentation_.

.. _documentation: http://resolwe-bio.readthedocs.org/

Install
=======

Prerequisites
-------------

Make sure you have Python_ (2.7 or 3.4+) installed on your system. If you don't
have it yet, follow `these instructions
<https://docs.python.org/3/using/index.html>`__.

Resolwe Bioinformatics requires PostgreSQL_ (9.4+). Many Linux distributions
already include the required version of PostgreSQL (e.g. Fedora 22+, Debian 8+,
Ubuntu 15.04+) and you can simply install it via distribution's package
manager. Otherwise, follow `these instructions
<https://wiki.postgresql.org/wiki/Detailed_installation_guides>`__.

.. _Python: https://www.python.org/
.. _PostgreSQL: http://www.postgresql.org/

Installing via pip
------------------

*Resolwe Bioinformatics is not available via PyPI yet.*

Installing from source
----------------------

This installation method will install all Resolwe Bioinformatics' dependencies
from PyPI_. Installing the ``psycopg2`` dependency will require having a C
compiler (e.g. GCC_) as well as Python and PostgreSQL development files
installed on the system.

.. note::

    The preffered way to install the C compiler and Python and PostgreSQL
    development files is to use your distribution's packages, if they exist.
    For example, on a Fedora/RHEL-based system, that would mean installing
    ``gcc``, ``python-devel``/``python3-devel`` and ``postgresql-devel``
    packages.

Download the `latest release of Resolwe Bioinformatics
<https://github.com/genialis/resolwe-bio/archive/master.tar.gz>`_ and extract
it.

Go to the directory with the extracted source code and install it::

    python setup.py install

.. _PyPi: https://pypi.python.org/
.. _GCC: https://gcc.gnu.org/

Contribute
==========

We welcome new contributors. To learn more, read Contributing_ section of the
documentation.

.. _Contributing: http://resolwe-bio.readthedocs.org/en/latest/contributing.html
