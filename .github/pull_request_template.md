## Checklist
<!--
Mark an `[x]` for completed items.
-->

* [ ] Update CHANGELOG.rst for each commit separately:
  * Pay attention to write entries under the "Unreleased" section.
  * Mark all breaking changes as "**BACKWARD INCOMPATIBLE:**" and put them
    before non-breaking changes.
  * If a commit modifies a feature listed under "Unreleased" section,
    it might be sufficient to modify the existing CHANGELOG entry from previous
    commit(s).
* [ ] Bump the process version:
  * **MAJOR version (first number)**: Backward incompatible changes (changes
    that break the api/interface). Examples: renaming the input/output, adding
    mandatory input, removing input/output...
  * **MINOR version (middle number)**: add functionality or changes in a
    backwards-compatible manner. Examples: add output field, add non-mandatory
    input parameter, use a different tool that produces same results...
  * **PATCH version (last number)**: changes/bug fixes that do not affect
    the api/interface. Examples: typo fix, change/add warning messages...
* [ ] All inputs are used in process.
* [ ] All output fields have a value assigned to them.

<!--
Read the following guidelines and remove them before opening the pull request.
-->

## Additional guidelines

### Processes ###

* Code that was once accepted should be kept like it is. It should be changed
  only when its functionality has changed. Small parts may be changed
  to fix bugs (by changing code only to make it nicer you might introduce a bug).
* Set the optimal number of cores that the tool can use (check if parallelization
  can be performed).
* Description of the tool should be self explanatory with links to the tool's
  homepage and publication.
* Write explanatory input labels and add short command line parameter names
  in square brackets (for example: Maximum length [``--maxlength``]).
* Where possible provide default values for the process input fields.
  If they are provided, then ``required = False`` is not needed.
* Add descriptive error and warning messages.
* Use ``Cmd()`` for each call in the process (call of a tool, a
  Bash command, a Python script or an R script) and handle return codes.
* Add ``self.progress()`` calls where appropriate.

### Commit and PR messages ###

* Each commit should be minimal (i.e. 1 change) and self-contained (including
  tests).
* Mark incomplete PRs with the [WIP] (Work In Progress) tag.
* Do not paste links to private repositories from the public one.

### Tests ###

* Add meaningful names of test files:
  * Genomes: ``genome_species.fa.gz``
  * Annotations: ``annotation_species.gtf.gz``
  * Adapters: ``adapters-name/type-source``.
* Specific files should be saved in a folder with the name of the tested tool.
  This folder should include separate subfolders for inputs and outputs.
* Keep test files' size small. If there is no other option, add the large
  test file to the ``files/large`` directory so it will be tracked with Git LFS.
* Pay attention to edge cases: what will happen if any data is
  missing, if the length of some list is 0, if a string is given where one
  expects an integer, if a file is empty, ...
* When testing workflows, try to verify workflow success by checking one
  of the last steps. Pick the one that is most likely to fail. Also test
  that all steps have 'OK' status.

<!--
Thanks!
-->
