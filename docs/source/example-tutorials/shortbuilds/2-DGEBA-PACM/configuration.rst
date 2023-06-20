.. _dgeba_configuration_file:

The Configuration File
======================

We previously described an example :ref:`configuration file for making a system of polystyrene <configuration_files>`, where all the directives were explained in detail.  Here we will only highlight the ``constituents`` directive.

We'll use a very small, stoichiometric system of only 200 DGEBA's and 100 PACM's.

.. code-block:: yaml

  constituents:
    DGE:
      count: 200
      symmetry_equivalent_atoms: [[C1,C2],[C3,C4],[O1,O2]]
      stereocenters: [C3]
    PAC:
      count: 100
      symmetry_equivalent_atoms: [[N1,N2],[C1,C2]]
      stereocenters: [C1]

Notice that the ``contituents`` directive is where we specify both symmetry-equivalent atoms and stereocenters in each monomer by referencing atom names.  DGEBA has three sets of symmetry-equivalent atoms and one symmetry-unique stereocenter, for example.  ``HTPolyNet`` uses information on symmetry-equivalent atoms to "expand" the set of reactions to include all possible combinations of reactants that symmetry allows.  In the case of a general diepoxy with both C1 and C2 atoms designated as reactive and a general diamine with both N1 and N2 atoms designated as reactive, ``HTPolyNet`` expands the *single* primary-to-seconary amine reaction that generates DGE~C1-N1~PAC, which is explicitly specified in the configuration file, into a set of four reactions with distinct products:

1. PAC~N1-C1~DGE
2. PAC~N1-C2~DGE
3. PAC~N2-C1~DGE
4. PAC~N2-C2~DGE

``HTPolyNet`` only does this because we stipulate that atoms C1 and C2 in DGE be treated as symmetric, and that atoms N1 and N2 on PACM be treated as symmetric.  Now, from the standpoint of chemical structure, these atoms *are* strictly symmetry-equivalent, since the molecules themselves are symmetric.  However, ``HTPolyNet`` does not recognize symmetry; the user can instead use symmetry-equivalence as a shortcut to automatically build all possible intermolecular interactions.  This guarantees that all possible reactions are templated correctly during a cure.

Note also that since PAC~N1-C1~DGE is the product of the primary-to-secondary amine reaction *and* a reactant in the secondary-to-tertiary reaction, ``HTPolyNet`` must also expand the secondary-to-tertiary amine reactions into a set of *eight* distinct reactions whose products are:

1. PAC~N1-C1~DGE~N1-C1~DGE
2. PAC~N1-C2~DGE~N1-C1~DGE
3. PAC~N2-C1~DGE~N2-C1~DGE
4. PAC~N2-C2~DGE~N2-C1~DGE
5. PAC~N1-C1~DGE~N1-C2~DGE
6. PAC~N1-C2~DGE~N1-C2~DGE
7. PAC~N2-C1~DGE~N2-C2~DGE
8. PAC~N2-C2~DGE~N2-C2~DGE

Note that each of these reactions indicates a product where *only* N1 *or* N2 of a PACM is reacted with *either* a C1 *or* C2 of two distinct DGEBAs.  Because the two nitrogens are so far apart on PACM, there is no need to parameterize an oligomer in which *both* N1 *and* N2 are simultaneously bound to DGEBA carbons.  So we see that during cure, ``HTPolyNet`` will search for atoms that satisfy one of *twelve* possible reactions, all because we specified *two* reactions among reactants that have several sets of symmetry-equivalent atoms *and* the product of the first reaction is a reactant in the second.

Stereocenters are used by ``HTPolyNet`` to generate enantiomers strictly for building initial liquid systems that are racemic by default.  That is, ``HTPolyNet`` assumes that any monomer structure provided by an input ``mol2`` or ``pdb`` file is a single enantiomer (if it has one or more stereocenters), and instead of copying only this structure to build the liquid, it instead generates *all* possible enantiomers by flipping stereocenters and then randomly selects among these to build a liquid.  ``HTPolyNet`` does brute force enumeration of enantiomers, meaning that a molecule with *N* stereocenters will have 2\ :sup:`N`\ 
enantiomers.

In the case of of DGEBA, two atoms are designated by stereocenters, and since they each also belong to distinct groups of symmetry-equivalent atoms, a DGEBA has four unique stereocenters.  This means 16 unique DGEBA diastereomers are used to build the liquid.  Two of those atoms are not strictly chiral in the activated form, since they each have two methyl ligands, but only one of those methyl ligands contains the reactive C1 (or by symmetry, C2) atom that will form a bond, at which point they *become* chiral.  Since we don't want polymerization to introduce handedness, we make sure we have a racemic mixture to begin with with respect to all *potentially* chiral atoms.

PACM has only two relevant stereocenters, which are the carbons in the cyclohexyl rings to which the amine nitrogens are attached, so four distinct PACM diastereomers are used in building the initial liquid.

Now we can turn to actually :ref:`running the build <dgeba_run>`.