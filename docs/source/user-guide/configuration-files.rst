Configuration Files
~~~~~~~~~~~~~~~~~~~

An HTPolyNet configuration file is a list of keyword:value pairs in ``YAML`` format.  The ``Library`` subpackage as a few example configuration files in the ``cfg`` directory.

The items in a configuration file break down into two major classes:

1. Items that specify run-time *parameters*.
2. Items that specify the system *chemistry*.

HTPolyNet parameters
''''''''''''''''''''

Below is a table of parameter keywords and descriptions for parameters that govern the overall execution of HTPolyNet.

General parameters
^^^^^^^^^^^^^^^^^^

===============================    ==============  =====================
Parameter                          Type            Description (default)
===============================    ==============  =====================
``Title``                          str             Descriptive title
``gmx``                            str             ``gmx`` command (``gmx``(default) or ``gmx_mpi``)
``gmx_options``                    quoted string   options to pass to all ``gmx`` calls ('-quiet -nobackup')
``gmx_mdrun``                      quoted string   ``mdrun`` command; defaults to ``gmx (options) mdrun``
``initial_density``                float           initial system density in kg/m3 (300.0)
``CURE_initial_search_radius``     float           initial capture radius in nm (0.5)
``CURE_radial_increment``          float           capture radius increment (0.25)
``desired_conversion``             float           desired fraction of possible crosslink bonds to form
``max_CURE_iterations``            int             maximum number of CURE iterations to run prior to reaching desired conversion
``late_treshhold``                 float           conversion above which reactions are all treated with probability 1.0
``charge_method``                  string          "gas" for Gasteiger; "bcc" for bcc; passed directly to antechamber
===============================    ==============  =====================

Other control parameters govern detailed aspects of the CURE algorithm.  These involve MD simulations performed immediately prior to and immediately after new bond addition to the topology, in order to relax those bonds.

Dragging parameters
^^^^^^^^^^^^^^^^^^^

Prior to introducing new bonds, one has the option of *"dragging"* atoms destined to be bonded to each other closer together in a series of dragging simulations.  The series is composed of stages, each of which involves three ``gmx mdrun`` calls: (1) a minimization; (2) an NVT relaxation; and (3) an NPT relaxation.  Soon-to-be-bonded atoms are connected by fictitious (type-6) harmonic bonds with equilibrium distances set at the current separation distances and relatively weak spring constants.  With each successive stage, the bond lengths are reduced and the spring constants increased until the desired separation distance and spring constant are achieved.  Dragging is optional.

===============================    ==============  =====================
Parameter                          Type            Description (default)
===============================    ==============  =====================
``max_drag_stages``                int             number of drag stages to perform
``drag_limit``                     float           minimum distance each separation should achieve (nm); 0.0 turns off dragging (0.0)
``drag_nvt_steps``                 int             number of MD steps for NVT relaxation during dragging (-2, signals ``gmx mdrun`` to use the value in the mdp file)
``drag_npt_steps``                 int             number of MD steps for NPT relaxation during dragging (-2, signals ``gmx mdrun`` to use the value in the mdp file)
===============================    ==============  =====================


Bond relaxation parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

*After* new bonds are formed and all other bonded interactions, atom types, and charges are mapped from each bond's appropriate template, a series of *bond relaxation* MD simulations are performed.  These are in all ways similar to the optional *dragging* simulations except for the fact that here, the actual chemical bond parameters are progressively brought to their correct values as specified in the GAFF.  Bond relaxation is *required* because most new bonds are much longer than they should be at equilibrium.

=================================    ==============  =====================
Parameter                            Type            Description (default)
=================================    ==============  =====================
``max_bond_relaxation_stages``       int             number of bond relaxation stages to perform
``max_bond_relaxation_increment``    float           maximum change in bond length parameters during relaxation (0.0; if set above 0.0, overrides ``max_bond_relaxation_stages``)
``relax_nvt_steps``                  int             number of MD steps for NVT relaxation 
``relax_npt_steps``                  int             number of MD steps for NPT relaxation 
=================================    ==============  =====================

Chemistry parameters
''''''''''''''''''''

The system chemistries and initial composition are specified by a set of inter-referential YAML entries.

Top-level chemistry parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

================================= ====          ===========
Parameter                         Type          Description
================================= ====          ===========
``initial_composition``           dict          keys: monomer names, values: numbers of molecules in system
``use_symmetry_equivalent_atoms`` list          monomers for which symmetry-equivalent atoms are used
``stereocenters``                 dict          keys: monomer names; values: list of names of chiral atoms
``reactions``                     list          reaction dicts, one per reaction
================================= ====          ===========

The ``initial_composition`` dictionary is how the initial extensive composition of the system is specified.  For example,

.. code-block:: yaml

    initial_composition: { MONA: 100, MONB: 200 }

specifies that the initial liquid should be composed of 100 ``MONA`` monomers and 200 ``MONB`` monomers.

The ``use_symmetry_equivalent_atoms`` list allows you to specify which molecules have internal symmetry that should be considered when enumerating all reactions.  Any molecule listed here is subjected to a short, 1000-K NVT MD simulation in vacuum from which an average interatomic distance matrix is calculated.  Using this distance matrix, symmetry-equivalent atoms are identified.  Each atom acquires the ``sea-idx`` attribute that contains its "symmetry index"; any two atoms with the same symmetry index are considered symmetry-equivalent.  For any atom name explicityly specified in a reaction or identified as a stereocenter in the input configuration, HTPolyNet will idenfity all symmetry-equivalent atoms and build out the full set of unique reactions and stereocenters from this information.  If you do not explicitly list a molecule in the ``use_symmetry_equivalent_atoms`` list, HTPolyNet assumes no such molecular symmetry exists, and only those reactions explicitly listed in the configuration file are considered.  So, if you have a diepoxy monomer but only reference one of the two reactive carbons in a reaction, no reactions involving the second carbon will take place UNLESS that second carbon is symmetry-equivalent to the one named AND the monomer appears in the ``use_symmetry_equivalent_atoms`` list.

The ``stereocenters`` dictionary is keyed by molecule name, and each entry is a list of one or more atom names which are declared as chiral centers.  By default, HTPolyNet will use this declaration to generate a racemic mixture of stereoisomers of any monomer in this list for the initial liquid simulation.  In an input mol2 file generated by, say, obabel, the atoms typically are not uniquely named.  **This means that you the user must edit any input mol2 file to give a unique name to any atom you will declare a stereocenter.**


Reaction dicts
^^^^^^^^^^^^^^

HTPolyNet expects one or more reaction templates to be defined in the configuration file.  A reaction is defined by the precise pairs of atoms that become new covalent bonds.  To precisely define each such pair, the reaction must also identify one or more reactant molecules.  Each reaction also names a single product molecule.  HTPolyNet will build oligomer templates using these reactions and then GAFF-parameterize them.  The parameterizations are used during CURE to re-type atoms and reset charges after each new bond is formed.

=================== =====  ===========
Keyword             Type   Description
=================== =====  ===========
``name``            str    descriptive name
``stage``           str    "cure" or "post-cure"
``probability``     float  probability that bond will form in one iteration if identified (1.0)
``reactants``       dict   keyword: reactant key, value: reactant molecule name
``product``         str    name of product molecule
``atoms``           dict   keyword: atom key, value: atom dict
``bonds``           list   list of bond dicts, one item per bond formed in reaction
=================== =====  ===========

Atom dicts and bond dicts are explained below.

Atom dicts
^^^^^^^^^^

An atom dict uniquely identifies an atom in a molecule by virtue of the molecule's name, the atom's residue number inside the molecule (starts at 1), the unique name of the atom in that residue, and its number of available sacrificial hydrogens (that is, the number of crosslink bonds it can participate in).  In an input mol2 file generated by, say, obabel, the atoms typically are not uniquely named.  **This means that you the user must edit any input mol2 file to give a unique name to any atom you will declare a reactive atom.**

=================== ====  ===========
Keyword             Type  Description
=================== ====  ===========
``reactant``        key   key to reactant in reactant dict to which this atom max_bond_relaxation_stages
``resid``           int   residue index in reactant molecular sequence to which this atom belongs (begins at 1)
``atom``            str   name of atom within that residue
``z``               int   number of available crosslink bonds for this atom
=================== ====  ===========

Bond dicts
^^^^^^^^^^

============= ======= ===========
Keyword       Type    Description
============= ======= ===========
``atoms``     list    list with the two atom keys the comprise the bond
``order``     float   bond order (currently not used; we let antechamber decide)
============= ======= ===========

An example reaction-dict is shown below:

.. code-block:: yaml

    name:     Primary-to-secondary-amine
    stage:    cure
    reactants: { 1: DFA, 2: FDE }
    product:   DFAFDE
    probability: 1.0
    atoms:
         A: { reactant: 1, resid: 1, atom: N1, z: 2 }
         B: { reactant: 2, resid: 1, atom: C1, z: 1 }
    bonds:
        - { atoms: [A, B], order: 1 }
    
This reaction-dict defines the reaction between an amine-containing molecule (DFA) and an epoxy-containing molecule (FDE) to produce an oligomer (DFAFDE).   It is instructive to read this YAML from bottom up.  There is one bond in the list of bond between atoms "A" and "B".  These are keys in the atoms dict right above.  "A" is the N1 atom of resid 1 of reactant 1, and "B" is the C1 atom of resid 1 of reactant 2.  The reactant keys are associated with molecule names in the reactants dict.  We see here that reactant 1 is DFA and reactant 2 is FDE.  

The values of the ``atom:`` keys in the atoms dict entries are atom names **as they appear in the mol2 file of the associated resid**.  In **this** case, both DFA and FDE are **monomers**: they are molecules with a single resid in their sequence. (Reactant and product molecules need not be monomers, but HTPolyNet must be able to trace any molecule back to its monomeric constituents by recursively traversing all reactions.) That is, this implies there is a unique atom named ``N1`` in the file ``DFA.mol2``, and a unique atom ``C1`` in the file ``FDE.mol2``.

If either reactant in any reaction is tagged as having more than one other atom symmetry-equivalent to any reactive atom, the original reaction is "symmetry-expanded" into as many new reactions as necessary to fully enumerate all possible reactions.  For instance, if N1 is symmetry equivalent to N2 in FDA, and C1 is symmetry equivalent to C2 in FDE (this is actually the case), then the "N1-C1" reaction above is automatically expanded into four symmetry-equivalent reactions: "N1-C1", "N2-C1", "N1-C2", and "N2-C2".  Each will generate a unique product: "DFAFDE", "DFAFDE-1", "DFAFDE-2", and "DFAFDE-3".  All products are necessary symmetry-equivalent, meaning all have the same set of atom-to-atomtype, bond-to-bondtype, etc., mappings, even though the particular atom indices are unique in each.
