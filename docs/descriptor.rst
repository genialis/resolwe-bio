==================
Descriptor schemas
==================

When working with the biological data, it is recommended (and often required) to
properly annotate samples. The annotation information attached to the samples
includes information about `organism`, `source`, `cell type`, `library preparation
protocols` and others.

The annotation fields associated with the samples or related sample files are
defined in the descriptor schemas. This tutorial describes the descriptor schemas
that are attached to the sample objects, raw sequencing reads and differential
expressions files.

Other available descriptor schemas can be explored at the Resolwe-bio GitHub_ page.
Customized descriptor schemas can be created using the `Resolwe SDK`_.


.. _GitHub: https://github.com/genialis/resolwe-bio/tree/master/resolwe_bio/descriptors
.. _Resolwe SDK: http://resdk.readthedocs.io/en/latest/index.html

Sample
======

When a new data object that represents a biological sample (i.e. fastq files,
bam files) is uploaded to the database, the unannotated sample ( presample) is
automatically created. When annotation is attached to the presample object, this
object is automatically converted to the annotated sample. To annotate the sample,
we need to define a descriptor schema that will be used for the annotation.
Together with the descriptor schema, we need to provide the annotations
(descriptors) that populate the annotation fields defined in the descriptor shema.
The details of this process are described in the `Resolwe SDK`_ documentation.

To annotate the sample in a GEO compliant way, we prepared the `sample`_
annotation schema. An example of the customized descriptor schema is also
`available`_.

.. _sample: https://github.com/genialis/resolwe-bio/blob/master/resolwe_bio/descriptors/sample_geo.yml
.. _available: https://github.com/genialis/resolwe-bio/blob/master/resolwe_bio/descriptors/sample_detailed.yml


Reads
=====

To annotate raw sequencing reads we have prepared two descriptor schemas: `reads`_
and `reads_detailed`_.

.. _reads: https://github.com/genialis/resolwe-bio/blob/master/resolwe_bio/descriptors/reads.yml
.. _reads_detailed: https://github.com/genialis/resolwe-bio/blob/master/resolwe_bio/descriptors/sample_detailed.yml



Differential expression
========================

To define the default thresholds for ``p-value``, ``log fold change (FC)``
and to describe which samples are used as cases and which as controls in
the calculation of differential expression we have prepared `diffexp`_
descriptor schema.

.. _diffexp: https://github.com/genialis/resolwe-bio/blob/master/resolwe_bio/descriptors/diffexp.yml