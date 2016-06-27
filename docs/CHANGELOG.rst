##########
Change Log
##########

All notable changes to this project are documented in this file.


==========
Unreleased
==========

Changed
-------
* Removed obsolete utilities superseded by resolwe-runtime-utils
* Use FDR values for Y axis in Limma volcano plot

Fixed
-----
* Sample descriptor schema updated
* Included all source files and supplementary package data in sdist
* Update gatk-joint processor

Added
-----
* "flow_collection: sample" to processes
* MACS14 processor
* Initial Tox configuration for running the tests
* Tox tests for ensuring high-quality Python packaging


================
1.0.0 2016-03-31
================

Changed
-------
* Renamed assertFileExist to assertFileExists
* Restructured processes folder hierarchy
* Removed re-require and hard-coded tools' paths

Fixed
-----
* Different line endings are correctly handled when opening gzipped files
* Fail gracefully if the field does not exist in assertFileExists
* Enabled processor tests (GO, Expression, Variant Calling)
* Enabled processor test (Upload reads with old Illumina QC encoding)
* Made Resolwe Bioinformatics work with Resolwe and Docker

Added
-----
* Import expressions from tranSMART
* Limma differential expression (tranSMART)
* VC filtering tool (Chemical mutagenesis)
* Additional analysis options to Abyss assembler
* API endpoint for Sample
* Initial documentation
