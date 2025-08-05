============
Contributing
============

Installing prerequisites
========================

Make sure you have Python_ 3.10+ installed on your system. If you don't have it
yet, follow `these instructions
<https://docs.python.org/3/using/index.html>`__.

.. _Python: https://www.python.org/

The pip_ tool will install all Resolwe Bioinformatics' dependencies from PyPI_.
Installing some (indirect) dependencies from PyPI_ will require having a C
compiler (e.g. GCC_) as well as Python development files installed on the
system.

.. note::

    The preferred way to install the C compiler and Python development files is
    to use your distribution's packages, if they exist.
    For example, on a Fedora/RHEL-based system, that would mean installing
    ``gcc`` and ``python3-devel`` packages.

Optional prerequisites
----------------------

Running Resolwe bio tests requires Docker_ to be installed on your system.

If you want to run or develop tests with large input or output files, then
install the `Git Large File Storage`_ extension.

.. _pip: https://pip.pypa.io/
.. _PyPi: https://pypi.python.org/
.. _GCC: https://gcc.gnu.org/
.. _Git Large File Storage: https://git-lfs.github.com/
.. _Docker: https://docs.docker.com/get-started/get-docker/

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

    pip install --pre -e .[docs,package,test]

.. note::

    We recommend using `pyvenv <http://docs.python.org/3/library/venv.html>`_
    to create an isolated Python environment for Resolwe Bioinformatics.

.. _Resolwe Bioinformatics' git repository: https://github.com/genialis/resolwe-bio

Running tests
=============

Manually
--------

Change directory to the ``tests`` Django project::

    cd tests

Run docker::

    docker compose up -d

.. note::
    On Mac or Windows, Docker might complain about non-mounted volumes.
    You can edit volumes in *Docker => Preferences => File Sharing*
    The following volumes need to be shared:

    - /private
    - /tmp
    - /var/folders


    ``/private`` is shared by default. When you attempt to add ``/var/folders``
    it might try to add ``/private/var/folders`` which will cause Docker complaining
    about overlapping volumes. Here's a workaround: Change ``/private`` to
    ``/var/folders`` and then add ``/private`` again.

Before running the tests, prepare the database::

    ./manage.py migrate
    ./manage.py createsuperuser

To run the tests, use::

    ./manage.py test resolwe_bio --parallel 2

.. note::

    If you don't specify the number of parallel test processes (i.e. you just
    use ``--parallel``), Django will run one test process per each core
    available on the machine.

.. warning::

    If you run Docker in a virtual machine (i.e. if you use MacOS or Windows)
    rather that directly on your machine, the virtual machine can become
    totally unresponsive if you set the number of parallel test processes too
    high. We recommend using at most ``--parallel 2`` in such cases.

To run a specific test, use::

    ./manage.py test resolwe_bio.tests.<module-name>.<class-name>.<method-name>

For example, to run the ``test_macs14`` test of the
``ChipSeqProcessorTestCase`` class in the ``test_chipseq`` module, use::

    ./manage.py test resolwe_bio.tests.processes.test_chipseq.ChipSeqProcessorTestCase.test_macs14

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

.. note::

    To control the number of test processes `Django will run in parallel`_, set
    the ``DJANGO_TEST_PROCESSES`` environment variable.

Since running tests for all processes may take a long time, there is an option
to run partial tests based on what files have been changed between HEAD and a
specific commit (e.g. master). The Tox environments that run partial tests have
the ``-partial`` suffix, e.g.::

    tox -e py312-partial

To configure the commit against which the changes are compared you should set
the ``RESOLWE_TEST_ONLY_CHANGES_TO`` environmental variable (it is set to master
by default).

.. _Tox: http://tox.testrun.org/
.. _Django will run in parallel:
    https://docs.djangoproject.com/en/1.10/ref/django-admin/#cmdoption-test-parallel

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

.. note::

    To build the documentation, you must use Python 3 (Python 2 is not
    supported).

Preparing release
=================

Follow `Resolwe's documentation on preparing a release`_. Resolwe
Bioinformatics code is automatically released to the PyPI when tagged.

.. _Resolwe's documentation on preparing a release:
  http://resolwe.readthedocs.io/en/latest/contributing.html#preparing-release
