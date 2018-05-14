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
* [ ] All output fields have their ``re-save`` calls.

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
* Description of tool should be self explanatory with links to tool's homepage
  and publication.
* Write explanatory input labels and add short command line parameter names
  in square brackets (for example: Maximum length [``--maxlength``]).
* Where possible provide default values for the process input fields.
  If they are provided, then ``required: false`` is not needed.
* Add descriptive error and warning messages.
* Place a ``re-checkrc`` call after each process step (call of a tool, a
  Bash command, a Python script or an R script).
* Add ``re-progress`` calls where appropriate.

### Commit and PR messages ###

* Each commit should be minimal (i.e. 1 change) and self-contained (including
  tests).
* Changes to Docker images need to be accepted in a separate pull request.
  Before the new images can be used for automated testing of your PR, they
  must be manually tagged and pushed to Docker Hub.
* Mark incomplete PRs with the [WIP] (Work In Progress) tag.
* Do not paste links to private repositories from the public one.

### Tests ###

* Add meaningful names of test files with spaces:
  * Genomes: ``genome species.fa.gz``
  * Annotations: ``annotation species.gtf.gz``
  * Adapters: ``adapters-name/type-source``.
* Output and specific input files should include the name of the tested tool.
* Keep test files' size small. If there is no other option, add the large
  test file to the ``files/large`` directory so it will be tracked with Git LFS.
* Pay attention to edge cases: what will happen if any data is
  missing, if length of some list is 0, if string is given where one expects
  integer, if file is empty, ...
* When testing workflows, try to verify workflow success by checking one
  of the last steps. Pick the one that is most likely to fail. Also test
  that all steps have 'OK' status.

<!--
Thanks!
-->
