Symmetry equivalent atoms
~~~~~~~~~~~~~~~~~~~~~~~~~

HTPolyNet allows for the use of symmetry-equivalent atoms in monomer structures.  Individual atoms are referenced in the configuration file as either reactive or a stereocenter, with reference to a particular molecule/residue.  If symmetry-equivalence is enabled for a molecule via the `use_symmetry_equivalent_atoms` option, then all atoms that are symmetry-equivalent to a specified atom are also considered.  Symmetry-equivalence is encoded in an atom coordinate attribute `sea-idx`, which contains an integer designation of the cluster of symmetry-equivalent atoms to which the attribute owner belongs.  All atoms with the same `sea-idx` (within one residue!) are considered symmetry-equivalent.  `sea-idx` attributes are stored in `*.sea` files for parameterized molecules.

What symmetry-equivalence does 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Sterocenters:  If an atom is designated a stereocenter, all atoms symmetry-equivalent to it are also designated stereocenters.  
2. Reactions:  Suppose a reaction specifies an atom `A` in reactant 1 and `X` in reactant 2.  Now suppose `A` is symmetry-equivalent to `B` in reactant 1, and `Y` is symmetry-equivalent to `X` in reactant 2.  So, instead of there only being a single reaction between 1 and 2 (where `A` bonds to `X`), there are now *four*: `A-X`, `A-Y`, `B-X`, and `B-Y`.  This results in four distinct, but symmetry-equivalent products.  When symmetry-equivalence is enabled, one need only specify the one reaction (`A-X`), and HTPolyNet will automatically generate the other three.

How symmetry-equivalence is determined
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For any molecule for which a `*.sea` file is not found, HTPolyNet will perform an "SEA" analysis.  This consists of performing a short, hot MD simulation of the molecule in vacuum and extracting from it the mean interatomic distance matrix.  In principle, two atoms are symmetry-equivalent if their distance-matrix columns comprise an identical set of distances.  In practice, for finite sampling from an MD simulation, two sets of distances are considered identical if, after putting the distances in order, the root-mean-squared distance between the columns is below some threshhold value.  Empirically, we have noticed that for small monomers of less than 60 or so atoms, a threshhold of 0.13 nm is sufficient to correctly identify the most topologically distance symmetry-equivalent heavy atoms, when a simulation is performed at 1000 K for 50,000 time-steps of 0.002 fs and sampling every 500 time-steps.

If you prefer not to automatically determine symmetry-equivalence this way (which is, we admit, not foolproof), then you can generate your own `*.sea` file for your monomer.