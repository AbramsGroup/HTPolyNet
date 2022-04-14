## Molecular structure data files for monomers and oligomers

* `monomers`  Files for monomers; must minimally include any user-required `mol2` files if these
are not provided by the user in the runtime directory.
* `oligomers` Files for oligomers, automatically generated when necessary.

Several file types are represented here

1. `mol2`  Sybyl MOL2 format files.  There are two types:
   * `XYZ.mol2` is interpreted as user input describing the atomic structure and topology (bonding information) of the `XYZ` monomer.  Some common monomer inputs are provided.  These files must conform well enough to Sybyl format to permit `antechamber` to read them.  
   * `XYZ-p.mol2` is the output of the parameterization process for monomer `XYZ`.  The atom names and types in this file are used in the simulations.
2. `gro`  Gromacs coordinate files generated automatically.  Usually in the form `XYZ-p.gro`.
3. `top`/`itp`  Gromacs topology files generated automatically.
4. `sea` Our custom symmetry-equivalent atom data; congruent to `gro` files, but not included in order to preserve the `gro` file format.  Symmetry-equivalent atoms for new monomers are detected automatically (if user desires) to reduce the number of oligomers that must be generated and parameterized.  `sea` files are only needed for monomers.

Additionally, if an end-capping chemistry is specified for any monomer at run-time, the files `XYZ-capped.mol2/gro/top/itp/sea` are also automaticaly generated and put in the `monomers` directory.
 