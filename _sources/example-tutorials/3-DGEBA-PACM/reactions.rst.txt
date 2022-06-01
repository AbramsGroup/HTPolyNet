.. _dgeba_reaction_dictionaries:

Reactions
=========

Now that the input structures for DGEBA and PACM have been generated, we are ready to generate the reaction dictionaries that describe the chemistry.  These will be stored in a YAML-format file called whatever you like, but it should reside in the top-level directory of the run (``my_dgeba_pacm_build``).

In a system of diepoxies and diamines, the two distinct crosslinking reactions are 

1. Reaction of an epoxy with a primary amine to generate a secondary amine; and
2. Reaction of an epoxy with a secondary amine to generate a tertiary amine.

In the case of PACM and DGEBA, these can be represented as

.. image:: r1.png

.. image:: r2.png

The first reaction generates a product we call "PACDGE", and the second reaction uses PACDGE as a reactant, along with DGE, to make "PACDGE2". Note that both result in deletion of two "sacrificial hydrogens". PACDGE serves as a parameterization template for new primary amines, and PACDGE2 serves as a parameterization template for new secondary amines.  That is, upon formation of a new secondary amine, the system acquires new bonded interactions, atoms types, and partial charges, based on those in the template that are contiguous with the new bond and new bonded interactions.

To allow for specification of general reactions, HTPolyNet parses a standard YAML syntax called a "reaction dictionary".  Here is the reaction dictionary for the first reaction:

.. code-block:: yaml

    name:     'Primary-to-secondary-amine'
    stage:     cure
    reactants: { 1: PAC, 2: DGE }
    product:   PACDGE
    probability: 1.0
    atoms: {
        A: {reactant: 1, resid: 1, atom: N1, z: 2}
        B: {reactant: 2, resid: 1, atom: C1, z: 1}
    }
    bonds:
        - {atoms: [A, B], order: 1}


Like all YAML, these keyword:value pairs.  ``name`` is just a label; you can use any meaningful string here.  ``stage`` is either ``cure`` or ``post-cure``; we'll see an example of a post-cure reaction below.  ``reactants`` labels the reactants with internal keys (here, ``1`` and ``2``), and the values ``PAC`` and ``DGE`` imply the existence of ``PAC.mol2`` and ``DGE.mol2``.  ``product`` names the product molecule.  ``probability`` is a number between 0 or 1 that governs the likelihood of the bond forming in any given iteration; that is, each bond matching this template that has passed all other single-bond and global filters is only allowed to form if its intrinsic probability is greater than a random number drawn between 0 and 1.  A probability of 1.0 implies matching bonds will always form.  ``atoms`` labels all atoms participating in the reaction, and we generally expect at least two atoms, one from each reactant.  Here, we label the N1 atom of PAC as atom ``A`` and the C1 atom of DGE as ``B``.  (We'll refer to these labels in just a minute.)  Notice that we reference PAC as "residue 1 of reactant 1", and we insist that the N1 atom have **two** available crosslink sites (``z`` is 2), and we reference DGE as "residue 1 of reactant 2", with a ``z`` of 1.  Finally, ``bonds`` is a list of bond dicts.  Each bond dict has two items: ``atoms`` is the list of the two atoms (referenced according to their keys in the ``atoms`` dict) that form the bond, and ``order`` is the bond order (this is currently unused).

Now, let's consider the reaction dictionary for the conversion of a secondary amine to a tertiary amine:

.. code-block:: yaml

    name:     'Secondary-to-tertiary-amine'
    stage:     cure
    reactants: { 1: PACDGE, 2: DGE }
    product:   PACDGE2
    probability: 0.5
    atoms: {
        A: {reactant: 1, resid: 1, atom: N1, z: 1}
        B: {reactant: 2, resid: 1, atom: C1, z: 1}
    }
    bonds:
        - {atoms: [A, B], order: 1}

Notice carefully the new reactant 1 here is PACDGE, and atom A is the N1 of residue 1 of reactant 1 with a ``z`` of 1.  This reaction generates the product PACDGE2.  We also set its probability to 0.5, which is approximately the same as saying these reactions intrisically happen with half the frequency of the secondary-amine-formation reactions.  This allows a way to specify the "relative reactivity" of the two types of reactions. (Typically, secondary-to-tertiary amine reactions are slower than analogous primary-to-secondary amine reactions.)

Finally, post-cure reactions that get rid of "extra" sacrificial hydrogens can be defined.  Consider:

.. image:: r3.png

This reaction removes the two sacrificial H's used to "open" the oxirane and reforms the oxirane.  The YAML syntax for this reaction is

.. code-block:: yaml

    name:     'Oxirane-formation'
    stage:     post-cure
    reactants: { 1: DGE }
    product:   DGEC
    probability: 1.0
    atoms: {
        A: {reactant: 1, resid: 1, atom: C1, z: 1}
        B: {reactant: 1, resid: 1, atom: O1, z: 1}
    }
    bonds:
        - {atoms: [A, B], order: 1}

With these reaction dictionaries defined, we are ready to create the :ref:`configuration file <dgeba_configuration_file>`.