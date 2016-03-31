======================
Resolwe Bioinformatics
======================

|docs|

.. |docs| image:: https://readthedocs.org/projects/resolwe-bio/badge/?version=latest
    :target: http://resolwe-bio.readthedocs.org/
    :alt: Documentation Status

Bioinformatics pipelines for the Resolwe_ dataflow package for `Django
framework`_.

.. _Resolwe: https://github.com/genialis/resolwe
.. _Django framework: https://www.djangoproject.com/


Docs & Help
===========

Read about getting started and how to write `processes` in the documentation_.

To chat with developers or ask for help, join us on Slack_.

.. _documentation: http://resolwe-bio.readthedocs.org/
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

.. _Contributing: http://resolwe-bio.readthedocs.org/en/latest/contributing.html
