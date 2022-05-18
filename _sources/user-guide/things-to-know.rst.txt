Important Things to Know
~~~~~~~~~~~~~~~~~~~~~~~~

1.  HTPolyNet uses the concept of "reactions" to describe crosslinking chemistry and any post-cure chemistry.  Reactions are specified in the input configuration file.
2.  HTPolyNet's "molecules" are the smallest chemical units required to represent (i) individual monomers, and (ii) crosslink bonds involving two or more monomers.  They are used for deriving atom types and parameters.  Molecules act as reactants or products in reactions, and initial composition is specified using molecule names.  Input molecular structures of monomers are required.
3. HTPolyNet uses the idea of "sacrificial hydrogens".  This means that it expects a monomer structure to represent a state in which it must lose an H to make an intermonomer bond.
4. HTPolyNet's reactions require specification of "reactive atoms" as those that participate in intermonomer bonds.  These atoms must have unique atom names in any input file.

