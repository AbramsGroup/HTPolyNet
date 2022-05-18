## Library package for HTPolyNet

This is a set of resources used by `HTPolyNet`, including example input configurations, Gromacs parameter files (`mdp`), and data files containing structural/topological information on monomers.

The `HTPolyNet/projectfilesystem` module manages this resource via its `Library` and `ProjectFileSystem` classes.  The `Library` resource is organized according to file extension and molecule class:

The file `XYZ.ext` resides in a subdirectory determined like this:

* If `ext` is `yaml`, `json`, or `cfg`,  then `XYZ.ext` is found in `cfg/`.
* If `ext` is `mdp`, the `XYZ.ext` is found in `mdp/`.
* If `ext` is one of `mol2`, `gro`, `top`, `itp`, or `sea`, then `XYZ.ext` is found *either* in  `molecules/monomers/` or in `molecules/oligomers/`.  `molecules/monomers/` has search preference.  Generally, filename prefixes for automatically generated oligomer file contain at least one `@` character.  (See documentation for the `oligomerize()` function in the `molecule` module.)


`Library` uses three main methods:

* `fetch()` by default copies to the current runtime directory, controlled by `HTPolyNet`.

* `store()` is the opposite of `fetch()`, and is used to store any automatically generated files for potential later use.  It shouldn't be necessary to parameterize the same molecule over and over.

* `exist()` simply checks to see that the file exists in its expected place in the library.