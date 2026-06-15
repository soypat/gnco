# _import-examples
This folder contains examples that need external imports to work.
In Go when you add a underscore prefix to a directory or file name it is
excluded from being compiled with the rest of the program.

This is done to avoid adding dependencies to gnco's go.mod file.
- Trivial Auditing
    - No dependency analysis needed for vulnerability scanning
- Ensures nothing outside gnco gets compiled maintaining binary sizes small
- Eliminate a whole set of attack vectors for gnco importers

