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

Optional prerequisites
----------------------

If you want to run or develop tests with large input or output files, then
install the `Git Large File Storage`_ extension.

.. _pip: https://pip.pypa.io/
.. _PyPi: https://pypi.python.org/
.. _GCC: https://gcc.gnu.org/
.. _Git Large File Storage: https://git-lfs.github.com/

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

Manually
--------

Change directory to the Django test project::

    cd tests

To run the tests, use::

    ./manage.py test resolwe_bio

To run a specific test, use::

    ./manage.py test resolwe_bio.tests.<module-name>.<class-name>.<method-name>

For example, to run the ``test_macs14`` test of the
``ChipSeqProcessorTestCase`` class in the ``test_chipseq`` module, use::

    ./manage.py test resolwe_bio.tests.test_chipseq.ChipSeqProcessorTestCase.test_macs14

Using Tox
---------

To run the tests with Tox_, use::

    tox

To re-create the virtual environment before running the tests, use::

    tox -r

To only run the tests with a specific Python version, use::

    tox -e py<python-version>

For example, to only run the tests with Python 3.5, use ::

    tox -e py35

.. note::

    To see the list of available Python versions, see ``tox.ini``.

.. _Tox: http://tox.testrun.org/

Running tests skipped on Docker
-------------------------------

To run the tests that are skipped on Docker due to failures and errors, set the
``RESOLWEBIO_TESTS_SKIP_DOCKER_FAILURES`` environment variable to ``no``.

For example, to run the skipped tests during a single test run, use::

    RESOLWEBIO_TESTS_SKIP_DOCKER_FAILURES=no ./manage.py test resolwe_bio

To run the skipped tests for the whole terminal session, execute::

    export RESOLWEBIO_TESTS_SKIP_DOCKER_FAILURES=no

and then run the tests as usual.

Running tests with large files
------------------------------

To run the tests with large input or output files, ensure you have the
`Git Large File Storage`_ extension installed and run the tests as usual.

Adding tests with large files
-----------------------------

If a test file is larger than 1 MiB, then put it in the
``resolwe_bio/tests/files/large/`` directory. Git Large File Storage
(LFS) extension will automatically pick it up and treat it appropriately.

To ensure contributors without Git LFS or users using the source distribution
can smoothly run the tests, decorate the tests using large files with the
following::

    @skipUnlessLargeFiles(<large-file1>, <large-file2>, ...)

where ``<large-file1>``, ``<large-file2>``, ... represent the names of large
files used inside a particular test.

The decorator will ensure the test is skipped unless these files are present
and represent real large files (not just Git LFS pointers).

Building documentation
======================

.. code-block:: none

    python setup.py build_sphinx

Preparing release
=================

Follow `Resolwe's documentation on preparing a release`_.

.. _Resolwe's documentation on preparing a release:
  http://resolwe.readthedocs.io/en/latest/contributing.html#preparing-release
