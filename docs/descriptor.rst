=========================================
Sample annotations and Descriptor schemas
=========================================

When working with the biological data, it is recommended (and often required) to
properly annotate samples and sample data.

Sample metadata / sample annotations are managed through the `Resolwe SDK`_. You can find
the tutorial on how to manage sample annotation is the `Annotate Samples`_ section of the
Resolwe SDK documentation.

The annotation fields associated with the sample data are defined in the descriptor schemas.
This tutorial describes the descriptor schemas that are attached to the raw sequencing reads
and differential expressions files. Other available descriptor schemas can be explored at the
Resolwe-bio GitHub_ page.


.. _GitHub: https://github.com/genialis/resolwe-bio/tree/master/resolwe_bio/descriptors
.. _Resolwe SDK: http://resdk.readthedocs.io/en/latest/index.html
.. _Annotate Samples: https://resdk.readthedocs.io/en/stable/tutorial-create.html#annotate-samples


Reads
=====

To annotate raw sequencing reads we have prepared corresponding `reads`_ descriptor schema.

.. _reads: https://github.com/genialis/resolwe-bio/blob/master/resolwe_bio/descriptors/reads.yml



Differential expression
========================

To define the default thresholds for ``p-value``, ``log fold change (FC)``
and to describe which samples are used as cases and which as controls in
the calculation of differential expression we have prepared `diffexp`_
descriptor schema.

.. _diffexp: https://github.com/genialis/resolwe-bio/blob/master/resolwe_bio/descriptors/diffexp.yml
