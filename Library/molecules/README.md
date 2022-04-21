## Molecular structure data files for monomers and oligomers

Several file types are represented here

1. `mol2`  Sybyl MOL2 format files.  There are two types:
   * `XYZ.mol2` is interpreted as _user input_ describing the atomic structure and topology (bonding information) of the `XYZ` molecule.  Some MOL2 inputs are provided for common monomers.  These files must conform well enough to Sybyl format to permit `antechamber` to read them.
   * `XYZ-p.mol2` is the output of the parameterization process for molecule `XYZ`.  The atom names and types in this file are used in the simulations.
2. `gro`  Gromacs coordinate files generated automatically.  Usually in the form `XYZ-p.gro`.
3. `top`/`itp`  Gromacs topology files generated automatically.    Usually in the form `XYZ-p.top/itp`.
4. `sea` Our custom symmetry-equivalent atom data; congruent to `gro` files, but not included inside `gro` files in order to preserve the `gro` file format.  Symmetry-equivalent atoms in any molecule are detected automatically if molecule names are found in the `use_sea` list of the configuration.

 