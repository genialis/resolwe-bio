=========
Type tree
=========

Process types are listed alphabetically. Next to each type is
a list of processes of that type. Types are hierarchical, with
levels of hierarchy separated by colon “:”. The hierarchy defines
what is accepted on inputs. For instance, `Expression (Cuffnorm)`_
process' input is ``data:alignment:bam``. This means it also
accepts all subtypes (*e.g.,* ``data:alignment:bam:bwasw``,
``data:alignment:bam:bowtie1`` and ``data:alignment:bam:tophat``).
We encourage the use of existing types in custom processes.

.. _Expression (Cuffnorm): catalog-definitions.html#process-upload-expression-cuffnorm

.. autoprocesstype::
