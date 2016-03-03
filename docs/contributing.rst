============
Contributing
============

Installing prerequisites
========================

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

The pip_ tool will install all Resolwe Bioinformatics' dependencies from PyPI_.
Installing the ``psycopg2`` dependency will require having a C compiler (e.g.
GCC_) as well as Python and PostgreSQL development files installed on the
system.

.. note::

    The preferred way to install the C compiler and Python and PostgreSQL
    development files is to use your distribution's packages, if they exist.
    For example, on a Fedora/RHEL-based system, that would mean installing
    ``gcc``, ``python-devel``/``python3-devel`` and ``postgresql-devel``
    packages.

.. _pip: https://pip.pypa.io/
.. _PyPi: https://pypi.python.org/
.. _GCC: https://gcc.gnu.org/

Preparing environment
=====================

`Fork <https://help.github.com/articles/fork-a-repo>`__ the main
`Resolwe Bioinformatics' git repository`_.

If you don't have Git installed on your system, follow `these
instructions <http://git-scm.com/book/en/v2/Getting-Started-Installing-Git>`__.

Clone your fork (replace ``<username>`` with your GitHub account name) and
change directory::

    git clone https://github.com/<username>/resolwe-bio.git
    cd resolwe-bio

Prepare Resolwe Bioinformatics for development::

    pip install -e .[docs,package,test]

.. note::

    We recommend using `virtualenv <https://virtualenv.pypa.io/>`_ (on
    Python 2.7) or `pyvenv <http://docs.python.org/3/library/venv.html>`_ (on
    Python 3.4+) to create an isolated Python environment for Resolwe.

.. _Resolwe Bioinformatics' git repository: https://github.com/genialis/resolwe-bio

Preparing database
==================

Add a postgres user::

    createuser -s -r postgres

Running tests
=============

To run the tests, use::

    cd tests
    ./manage.py test resolwe_bio

To run the tests with Tox_, use::

    tox -r

.. _Tox: http://tox.testrun.org/

Building documentation
======================

.. code-block:: none

    python setup.py build_sphinx
