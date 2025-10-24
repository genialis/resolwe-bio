.. index:: Writing processes

=================
Writing processes
=================

Resolwe bioinformatics (resolwe-bio) contains a collection of bioinformatics
processes written for the `Resolwe`_ dataflow engine. `Resolwe processes`_
are defined by their input and output fields, process requirements
(e.g., execution environment, compute requirements), and the processing logic
that maps the inputs into outputs. Resolwe processes are written in Python and
can be joined into bioinformatics pipelines (workflows). An example
`Python process template`_ is available in the Resolwe documentation.

.. _Resolwe: https://resolwe.readthedocs.io/en/latest/
.. _Resolwe processes: https://resolwe.readthedocs.io/en/latest/proc.html
.. _Python process template: https://github.com/genialis/resolwe/blob/master/docs/example/example/processes/template_py_process.py

Example process
===============

Here we list an example process that imports paired-end FASTQ files to a Resolwe server.
The process accepts both mate1 and mate2 input FASTQ files and an optional
sample species annotation as input. The process stores the imported FASTQ files
as output and, if provided, also stores the `species` annotation in the sample
metadata.

.. literalinclude:: example/processes/import_fastq.py

.. note::

    The example process above is a simplified version of the
    `upload-fastq-paired` process included with Resolwe Bioinformatics. The actual
    process contains additional input and output fields, as well as additional
    logic in the main process to handle different input scenarios, file validation
    logic and QC steps.

The Resolwe process class consists of several fields, which are described in
detail in the `Resolwe processes`_ documentation. Following the example above, the
fields used in the process are:

- ``slug = "upload-fastq-paired-docs"``: a unique identifier of the process

- ``name = "Upload paired-end FASTQ files"``: a human-readable name of the process

- ``process_type = "data:reads:fastq:paired"``: a string that classifies the process

- ``version = "1.0.0"``: a version of the process. The version should be updated
  when the process code or any of its dependencies change in a way that
  affects the output of the process. Only the highest version of a process is
  available for execution on the Resolwe server.

- ``category = "Import"``: a string that classifies the process (e.g., 'Import')

- ``data_name = '{{ mate1.file|default("?") }}'``: a string that defines how
  the data object created by the process will be named. The string can contain
  Jinja2 template expressions that refer to input fields of the process.
  In this example, the data object will be named after the name of the mate1 FASTQ file,
  or "?" if mate1 input field would not be defined.

- ``scheduling_class = SchedulingClass.BATCH``: a scheduling class that defines
  how the process will be scheduled for execution on the Resolwe server. Possible values are
  ``SchedulingClass.BATCH`` (default) and ``SchedulingClass.INTERACTIVE``.
  Interactive processes are processed in a dedicated processing queue.

- ``persistence = Persistence.RAW``: a persistence class that defines how
  the data object created by the process will be treated in terms of data retention
  and reproducibility. Possible values are ``Persistence.RAW`` (default),
  ``Persistence.CACHED``, and ``Persistence.TMP``. ``Persistence.RAW`` type is
  to be used by the processes that import data files. ``Persistence.CACHED`` and
  ``Persistence.TMP`` data objects should be used for idempotent processing jobs,
  i.e. re-running the process with the same input should produce the same output.
  ``Persistence.CACHED`` data objects are stored permanently, while ``Persistence.TMP``
  data objects are temporary and can be deleted by the system when needed.

- ``entity_type = "sample"``: a field that defines the type of the entity
  to which the data object created by the process will be associated. In this case,
  the data object will be associated with a new ``sample`` object. When the Data object
  created by this process is used as an input to another process also marked with
  ``entity_type = "sample"``, the data object will be associated with the same
  sample object as the input data object.

- ``requirements = { ... }``: a dictionary that defines the process requirements
  (e.g., execution environment, compute requirements). The ``resources`` section
  defines the compute requirements of the process, including ``cores`` (integer), ``memory``
  (in MB, integer), and ``storage`` (in GB, integer). The ``docker`` section defines the
  execution environment of the process. Docker images commonly used by the Resolwe
  processes are available on `AWS ECR Public Gallery`_, and defined in the
  `resolwe-docker-images`_ repository. Custom Docker images can be built and
  pushed to a Docker registry, and then used by the processes.

``Input`` and ``Output`` classes of the Resolwe process contain the input and output
field definitions, respectively. The ``run`` method contains the processing logic
that maps the input fields to output fields. The processing logic is written in Python
and can use any Python packages and installed tools available in the execution environment
defined in the ``docker`` section of the process requirements.

The ``upload-fastq-paired-docs`` process uses the ``import_file()`` utility function
to import the input FASTQ files in compressed format. The returned file handles are
then evaluated and renamed to contain the ``fastq.gz`` suffix before they are assigned to
the output fields ``mate1`` and ``mate2``. If the optional ``species`` input field is
provided, the ``species`` annotation value is stored both to the dedicated data object output
field as well as to the sample annotation field ``species`` belonging to the ``general``
annotation field group.

.. _AWS ECR Public Gallery: https://gallery.ecr.aws/
.. _resolwe-docker-images: https://github.com/genialis/resolwe-docker-images

Tools
=====

Sometimes, it is very useful to write a custom script in R or other languages
to perform a certain task in process' algorithm. Custom scripts needed by
processes included with Resolwe Bioinformatics are located in the
`resolwe_bio/tools`_ directory.

.. note::

    A Resolwe's :class:`Flow Executor <resolwe.flow.executors.BaseFlowExecutor>`
    searches for tools in a Django application's ``tools`` directory or
    directories specified in the ``RESOLWE_CUSTOM_TOOLS_PATHS`` Django setting.

.. _resolwe_bio/tools: https://github.com/genialis/resolwe-bio/tree/master/resolwe_bio/tools

Writing workflows
=================

Workflows are collections of processes connected together to form a dataflow.
As processes, workflows are defined in Python syntax. We will create an example
workflow that produces the QC analysis of the FASTQ reads imported using the
``upload-fastq-paired-docs`` process described above. The workflow uses the
``fastqc-paired-end`` process defined in `fastqc.py`_ to perform the FastQC analysis.

.. _fastqc.py: https://github.com/genialis/resolwe-bio/tree/master/docs/example/processes/fastqc.py

.. literalinclude:: example/processes/qc_workflow.py

Writing and running tests
=========================

Processes and workflows included with Resolwe Bioinformatics are tested using unit tests.
These tests are located in the `resolwe_bio/tests`_ directory and are organized by process
and workflow types.

.. _resolwe_bio/tests: https://github.com/genialis/resolwe-bio/tree/master/resolwe_bio/tests

To run the tests, you need to have a working Resolwe-bio development environment with Docker
dependencies installed. Then, navigate to the ``resolwe-bio/tests`` directory and run:

.. code-block:: bash

    # Start the required docker containers
    docker compose up -d

.. note::

    By default, each test run triggers the download of all the required Docker images.
    To avoid this, and download only the images required by the tests you are running,
    you can set the ``RESOLWE_DOCKER_DONT_PULL`` environment variable to ``1`` in the
    shell environment, or directly in the ``tests/settings.py`` file.

A workflow test that runs the file upload process and triggers the example QC workflow
is shown below. The test uses the ``KBBioProcessTestCase`` class that provides
utilities for testing bioinformatics processes and workflows.

.. literalinclude:: example/tests/workflow_test.py

The test is triggered by calling the ``test`` management command on the specified test class.

.. code-block:: bash

    # Example test command. Adjust the path to the test file as needed.
    # This command runs the QC workflow test defined in the example workflow test file.
    # The input files used in the test should be placed in the
    # resolwe_bio/tests/files/ directory. The process scripts are expected to be
    # placed in the resolwe_bio/processes/ directory for the Resolwe server to successfully
    # locate and register them.
    ./manage.py test docs.example.tests.workflow_test.DocsProcessTestCase.test_qc_workflow

