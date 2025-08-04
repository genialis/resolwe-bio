.. index:: Writing processes

=================
Writing processes
=================

A short tutorial about writing bioinformatics pipelines (process is a
step in the pipeline) is contained in `Resolwe documentation`_. While the
Resolwe documentation describes the syntax of the YAML-based processes that are
slowly being deprecated, the Resolwe Bioinformatics processes are most often written
in Python. An example `Python process template`_ is available in the Resolwe documentation.

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

.. _Resolwe documentation: https://resolwe.readthedocs.io/en/latest/proc.html
.. _resolwe_bio/tools: https://github.com/genialis/resolwe-bio/tree/master/resolwe_bio/tools
.. _Python process template: https://github.com/genialis/resolwe/blob/master/docs/example/example/processes/template_py_process.py
