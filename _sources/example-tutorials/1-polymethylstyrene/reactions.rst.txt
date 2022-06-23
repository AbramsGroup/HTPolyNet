.. _pms_reaction_dictionaries:

Reactions
=========

In this system, the only *actual* reaction is between  a ``C1`` of one 4-methylstyrene monomer and the ``C2`` of another:

.. code-block:: yaml

 - {
      name:        'EMB1_1',
      stage:       cure,
      reactants:   {1: EMB, 2: EMB},
      product:     EMB1_1,
      probability: 1.0,
      atoms: {
        A: {reactant: 1, resid: 1, atom: C1, z: 1},
        B: {reactant: 2, resid: 1, atom: C2, z: 1}
      },
      bonds: [
        {atoms: [A, B], order: 1}
      ]
    }

Naturally, this generates a dimer:

.. image:: EMB1_1_labelled.png

Notice that atom ``C1`` of resid 1 has "attacked" ``C2`` of resid 2 to form the bond.

Now, consider the new bonded interactions that spring from the creation of this one bond.  Clearly, it is one side of a bunch of new angles, each of which is defined by one other atom on either resid 1 or 2.  It is also the central bond of many new proper dihedrals, each of which involves one neighbor atom each on resid 1 and another on resid 2.  Clearly, the dimer is sufficient to permit parameterization to yield these interactions.

However, that is not all.  The new bond is also one of the *outer* bonds of **many** dihedrals.  For example, the 1-C1--2-C2 bond is part of a torsion in which the 1-C2--1-C1 bond is the *central* bond, meaning each other atom on ``C2`` of resid 1 defines a unique proper dihedral.  In this dimer, all those atoms are hydrogens.  However, what if ``C2`` of resid 1 is attacked by ``C1`` of a third monomer?  Then *one* of these dihedrals is actually a ``C-C-C-C`` dihedral, not an ``H-C-C-C`` dihedral.  This means that, generally, for systems with monomers that react via double-bond opening, we need oligomeric templats that involve 3 and 4 monomers to have a full set of templates.

Naturally, a trimer and quadrimer suffice, so let's make those templates using zero-probability reactions:

.. code-block:: yaml

  - {
      name:        'EMB2_1',
      stage:       cure,
      reactants:   {1: EMB1_1, 2: EMB},
      product:     EMB2_1,
      probability: 0.0,
      atoms: {
        A: {reactant: 1, resid: 2, atom: C1, z: 1},
        B: {reactant: 2, resid: 1, atom: C2, z: 1}
      },
      bonds: [
        {atoms: [A, B], order: 1}
      ]
    }
  - {
      name:        'EMB1_2',
      stage:       cure,
      reactants:   {1: EMB, 2: EMB1_1},
      product:     EMB1_2,
      probability: 0.0,
      atoms: {
        A: {reactant: 1, resid: 1, atom: C1, z: 1},
        B: {reactant: 2, resid: 1, atom: C2, z: 1}
      },
      bonds: [
        {atoms: [A, B], order: 1}
      ]
    }
  - {
      name:        'EMB2_2',
      stage:       cure,
      reactants:   {1: EMB1_1, 2: EMB1_1},
      product:     EMB2_2,
      probability: 0.0,
      atoms: {
        A: {reactant: 1, resid: 2, atom: C1, z: 1},
        B: {reactant: 2, resid: 1, atom: C2, z: 1}
      },
      bonds: [
        {atoms: [A, B], order: 1}
      ]
    }

Note that this makes the trimer twice to account for whether the third monomer is a ``C1`` to ``C2`` attack, or vice versa.  Although the final molecules are the same, the particular reaction defines the template mapping, so it is safer (though maybe a bit redundant) to treat these two as "unique" templates.

The only time the quadrimer is needed is when a "chain"-"chain" reaction occurs such that dihedrals that must be templated comprise a set of atoms spread across four contiguous opened double-bonds.

Finally, we can choose to revert any fully unreacted monomers back to true double bonds:

.. code-block:: yaml

 - {
      name:         'EMBCC',
      stage:        post-cure,
      reactants:    {1: EMB},
      product:      EMBCC,
      probability:  1.0,
      atoms: {
        A: {reactant: 1, resid: 1, atom: C1, z: 1},
        B: {reactant: 1, resid: 1, atom: C2, z: 1}
      },
      bonds: [
        {atoms: [A, B], order: 2}
      ]
    }

The next thing we consider is the :ref:`configuration file <pms_configuration_file>` necessary to describe the crosslinking chemistry and determine the system build.