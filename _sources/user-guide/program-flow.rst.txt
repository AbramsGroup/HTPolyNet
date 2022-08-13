Program Flow
------------

.. figure:: flow1.png
   :scale: 80 %
   :alt: outer program flow

   ``htpolynet run`` workflow.

A basic depiction of the the workflow initiated by ``htpolynet run`` is shown in the figure above.  The "0"th step is generation of the molecular structure data for any monomeric reactants, in the form of either Sybyl ``mol2`` or RCSB ``pdb`` files.  ``HTPolyNet`` does **not** do this, we provide some general guidance :ref:`here <molecular_structure_inputs>` and some specific example cases in the :ref:`tutorials <example_tutorials>`.  Then, based on instructions in the :ref:`configuration file <configuration_files>`, ``HTPolyNet`` proceeds with setting up all reactions and oligomer templates.  Once these are generated, it then generates the full initial system topology in Gromacs format (that is, it generates a ``top`` file), and an initial set of coordinates (a ``gro`` file).  Since the initial coordinates are built a low density (typically), ``HTPolyNet`` then performs MD to "densify" the system, followed by any pre-cure equilibration the user would like.  Then the :ref:`CURE algorithm <cure_section>` takes over to generate intermolecular bonds and drive the polymerization.  Once it finishes, ``HTPolyNet`` conducts any post-cure equilibration the user likes, before saving the final ``top`` and ``gro`` file.  All of this work happens in the project subdirectory of the current directory in which `` run`` is invoked.

This workflow should make clear that the two required tasks of the user are:

1. :ref:`Generating monomer structure files; <molecular_structure_inputs>` and
2. :ref:`Creating an input configuration file. <configuration_files>`

.. _cure_section:

Connect-Update-Relax-Equilibrate (CURE) algorithm
-------------------------------------------------

.. figure:: CURE.png 
   :scale: 80 %
   :alt: picture of cure algorithm block flow diagram

   Block flow diagram of the CURE algorithm used in ``HTPolyNet``.

The algorithm used to create new bonds and polymerize a system is called the CURE algorithm, depicted above.  This is just a slightly modified version of a standard search-radius-type algorithm, first used by Li and Strahan to study EPON/DETDA thermosets (`Li and Strahan, Polymer 2010;51,6058 <https://doi.org/10.1016/j.polymer.2010.10.033>`_).  The CURE algorithm begins by executing a search for new bonds on a frozen system configuration.  Bonds are downselected through a series of filters to arrive at a final set of bonds to form.  If the distance between any pair of "bond-designate" atoms is greater than some threshold (the ``trigger_distance`` parameter in the :ref:`drag subdirective <cure.drag>` of the ``CURE`` directive of a configuration file), a series of MD simulations that slowly bring all to-be-bound atom closer together is performed.  Then the topology is updated, where ``HTPolyNet`` applies the charges, atom type, and bonded interaction templates from the oligomer template set to each bond.  After the update, a series of relaxation MD simulations bring all bonds to their equilibrium lengths.  Then a short NPT MD simulation equilibrates the overall density before initiating the next CURE iteration.  CURE iterations continue until (a) a desired conversion is reached, or (b) no new allowable bonds are identified.

.. _bondsearch_filters:

Bondsearch filters
------------------

A key task of ``HTPolyNet`` is identifying potential bonds between reactive atoms.  ``HTPolyNet`` organizes its "bondsearch" in any one CURE iteration along reactions.  Each reaction specified in the configuration file usually defines one bond by the identies of each reactant and the name of each atom in each reactant.  The bondsearch algorithm begins considering each reaction, and for each one, filtering the set of *possible* bonds using a series of rules, depicted below.

.. figure:: bond_filter.png 

   Bond filtering in one CURE iteration.

