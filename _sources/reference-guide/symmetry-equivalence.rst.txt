.. _symmetry_equivalence:

Symmetry equivalent atoms
~~~~~~~~~~~~~~~~~~~~~~~~~

HTPolyNet allows for the use of symmetry-equivalent atoms in monomer structures.  Individual atoms are referenced in the configuration file as either reactive or a stereocenter, with reference to a particular molecule/residue.  If symmetry-equivalence is enabled for a molecule via the ``symmetry_equivalent_atoms`` dictionary, then all atoms that are symmetry-equivalent to a specified atom are also considered.  An entry in the ``symmetry_equivalent_atoms`` is keyed by the molecule name, and the value is the list of "symmetry sets", each of which is a list of atom names.  One entry can have multiple symmetry sets, and each set must meaningfully have at least two atom names.  Symmetry-equivalence is encoded in an atom coordinate attribute `sea-idx`, which contains an integer designation of the cluster of symmetry-equivalent atoms to which the attribute owner belongs.  All atoms with the same `sea-idx` (within one residue!) are considered symmetry-equivalent.  `sea-idx` attributes are stored in `*.sea` files for parameterized molecules.

What symmetry-equivalence does 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Sterocenters:  If an atom is designated a stereocenter, all atoms symmetry-equivalent to it are also designated stereocenters.  
2. Reactions:  Suppose a reaction specifies an atom `A` in reactant 1 and `X` in reactant 2.  Now suppose `A` is symmetry-equivalent to `B` in reactant 1, and `Y` is symmetry-equivalent to `X` in reactant 2.  So, instead of there only being a single reaction between 1 and 2 (where `A` bonds to `X`), there are now *four*: `A-X`, `A-Y`, `B-X`, and `B-Y`.  This results in four distinct, but symmetry-equivalent products.  When symmetry-equivalence is enabled, one need only specify the one reaction (`A-X`), and HTPolyNet will automatically generate the other three.
