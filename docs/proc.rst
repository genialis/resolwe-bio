.. index:: Writing processes

=================
Writing processes
=================

A tutorial about writing bioinformatics pipelines (process is a
step in the pipeline) is in `Resolwe SDK for Python documentation`_.

Tools
=====

Frequently, it is very useful to write a custom script in Python or R to
perform a certain task in process' algorithm. For an example, see the
tutorial in `Resolwe SDK for Python documentation`_.

Custom scripts needed by processes included with Resolwe Bioinformatics are
located in the `resolwe_bio/tools`_ directory.

.. note::

    A Resolwe's :class:`Flow Executor <resolwe.flow.executors.BaseFlowExecutor>`
    searches for tools in a Django application's ``tools`` directory or
    directories specified in the ``RESOLWE_CUSTOM_TOOLS_PATHS`` Django setting.

.. _Resolwe SDK for Python documentation: http://resdk.readthedocs.io/en/latest/tutorial.html
.. _resolwe_bio/tools: https://github.com/genialis/resolwe-bio/tree/master/resolwe_bio/tools
