##########
Change Log
##########

All notable changes to this project are documented in this file.


==========
Unreleased
==========

Changed
-------
* Renamed assertFileExist to assertFileExists.
* Restructured processes folder hierarchy

Fixed
-----
* Different line endings are correctly handled when opening gzipped files.
* Fail gracefully if the field does not exist assertFileExists.
* Enabled processor tests (GO, Expression, Variant Calling)
* Enabled processor test (Upload reads with old Illumina QC encoding)

Added
-----
* Import expressions from tranSMART.
* Limma differential expression (tranSMART)
* Added VC filtering tool (Chemical mutagenesis)
