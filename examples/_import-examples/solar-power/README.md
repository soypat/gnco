# import-examples
This folder contains examples that need external imports to work. 

This is done to avoid adding dependencies to gnco's go.mod file.
- Trivial Auditing
    - No dependency analysis needed for vulnerability scanning
- Eliminate a whole set of attack vectors for gnco importers
- Ensures nothing outside gnco gets compiled maintaining binary sizes small
