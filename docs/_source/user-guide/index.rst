User Guide
==========

Overview
~~~~~~~~

Basic Usage
^^^^^^^^^^^

Installation of the HTPolyNet package gives access to the ``htpolynet`` command:

.. code-block:: console

    $ htpolynet -h 
    usage: htpolynet [-h] [-cfg CFG] [-lib LIB] [-log LOG] [-restart] [--force-parameterization] [--force-sea-calculation] [--force-checkin] [--loglevel LOGLEVEL] command

    positional arguments:
      command               command (info, run)

    optional arguments:
     -h, --help            show this help message and exit
     -cfg CFG              input config file (required)
     -lib LIB              local library, assumed flat (optional)
     -log LOG              log file (htpolynet_runtime.log)
     -restart              restart in latest proj dir
     --force-parameterization
                           force GAFF parameterization of any input mol2 structures (False)
     --force-sea-calculation
                           force calculation of symmetry-equivalent atoms in any input mol2 structures (False)
     --force-checkin       force check-in of any generated parameter files to the system library (False)
     --loglevel LOGLEVEL   Log level; info, debug (info)


The ``info`` command gives some basic information of your HTPolyNet installation:

.. code-block:: console

    $ htpolynet info
    This is some information on your installed version of HTPolyNet
    Library System Directories:
    /home/cfa/Git/HTPolyNet/Library/cfg
    /home/cfa/Git/HTPolyNet/Library/mdp
    /home/cfa/Git/HTPolyNet/Library/molecules
    /home/cfa/Git/HTPolyNet/Library/molecules/inputs
    /home/cfa/Git/HTPolyNet/Library/molecules/parameterized

    Commands available for HTPolyNet to use:
    antechamber (ver.   21.0) at /home/cfa/anaconda3/envs/py39/bin/antechamber     
        tleap (ver.   21.0) at /home/cfa/anaconda3/envs/py39/bin/tleap           
        parmchk2 (ver.   21.0) at /home/cfa/anaconda3/envs/py39/bin/parmchk2        
            gmx (ver. 2021.2) at /usr/local/gromacs/bin/gmx                        

The "Library" is a subpackage that provides some necessary and example inputs.  In particular, template ``mdp`` files for Gromacs are in the ``mdp`` directory.  In the ``molecules`` directory, the ``inputs`` subdirectory has a selection of ``mol2`` files for various monomers, while the ``parameterized`` subdirectory has ``top``, ``itp``, and ``sea`` files for molecules that have been previously parameterized during HTPolyNet development.

The Basic Algorithm
^^^^^^^^^^^^^^^^^^^

HTPolyNet produces Gromacs topology and coordinate files for amorphous, crosslinked polymer systems following the basic steps below:

1. If necessary, GAFF-parameterize any input monomers, and build any intermediate molecules dictated by the polymerization chemistry and GAFF-parameterize those.  Outputs of parameterizations are typically saved in a Library so that repeated parameterizations are not done.
2. Based on a specified composition, generate an initial simulation box, and equilibrate it to a liquid-like density.
3. Perform CURE (Connect-Update-Relax-Equilibrate) iterations to introduce intermonomer bonds according to the input chemistry until a desired conversion is met or an execution threshhold is reached.
4. Finalize by performing any requested post-cure chemistry.

Important Things to Know
^^^^^^^^^^^^^^^^^^^^^^^^

1.  HTPolyNet uses the concept of "reactions" to describe crosslinking chemistry and any post-cure chemistry.  Reactions are specified in the input configuration file.
2.  HTPolyNet's "molecules" are the smallest chemical units required to represent (i) individual monomers, and (ii) crosslink bonds involving two or more monomers.  They are used for deriving atom types and parameters.  Molecules act as reactants or products in reactions, and initial composition is specified using molecule names.  Input molecular structures of monomers are required.
3. HTPolyNet uses the idea of "sacrificial hydrogens".  This means that it expects a monomer structure to represent a state in which it must lose an H to make an intermonomer bond.
4. HTPolyNet's reactions require specification of "reactive atoms" as those that participate in intermonomer bonds.  These atoms must have unique atom names in any input file.

Inputs
~~~~~~

HTPolyNet requires minimally two types of inputs:

1. Molecular structures of any monomers as Tripos mol2 files; these are used directly as inputs to ``antechamber``.
2. An input configuration file in ``YAML`` format for controlling the ``htpolynet`` run.

Molecular structures
^^^^^^^^^^^^^^^^^^^^

HTPolyNet uses the distinct terms "monomer" and "molecule".  A monomer is akin to a "residue" for biomolecules; it has a unique name, and a "molecule" is made of a sequence of one or more monomers.  (Here, "sequence" only refers to the order in which monomers appear in the list of all atoms; not their topological sequence.)  HTPolyNet will build a system out of molecules you specify as input; typically, these would be monomers (i.e., molecules with a single monomer).

Regardless, any molecule you wish to use as a monomer must have a mol2 file ``<NAME>.mol2``, where ``<NAME>`` is replaced with the name of the monomer. This name is important; it is how the monomer is called forever inside HTPolyNet.

Sample mol2 files for a few monomers are provide in the Library subpackage: ``Library/molecules/inputs``.  You will likely want to create your own.  Most chemical structure drawing programs will output mol2 files.  OpenBabel can also generate them from SMILES strings; e.g., for butane:

.. code-block:: console

    $  echo "CCCC" | obabel -ismi -omol2 --gen3d --title BUTANE | sed s/"UNL1  "/"BUTANE"/ > BUTANE.mol2

So to generate mol2 files for your monomers, you have a few options.

The mol2 format HTPolyNet requires for a monomer is minimal and requires only three sections: ``@<TRIPOS>MOLECULE``, ``@<TRIPOS>ATOM`` and ``@<TRIPOS>BOND``.  For a molecule with more than one monomer, a ``@<TRIPOS>SUBSTRUCTURE`` section is also required.  Below is the mol2 file for butane we just created above::

    @<TRIPOS>MOLECULE
    BUTANE
    14 13 0 0 0
    SMALL
    GASTEIGER

    @<TRIPOS>ATOM
         1 C           0.9926    0.1046    0.0014 C.3     1  BUTANE     -0.0653
         2 C           2.5129    0.0970   -0.0099 C.3     1  BUTANE     -0.0562
         3 C           3.0598   -0.7349   -1.1682 C.3     1  BUTANE     -0.0562
         4 C           4.5801   -0.7416   -1.1802 C.3     1  BUTANE     -0.0653
         5 H           0.6228    0.7054    0.8382 H       1  BUTANE      0.0230
         6 H           0.5960   -0.9100    0.1102 H       1  BUTANE      0.0230
         7 H           0.5955    0.5315   -0.9250 H       1  BUTANE      0.0230
         8 H           2.8777   -0.3052    0.9424 H       1  BUTANE      0.0263
         9 H           2.8771    1.1281   -0.0871 H       1  BUTANE      0.0263
        10 H           2.6944   -0.3331   -2.1205 H       1  BUTANE      0.0263
        11 H           2.6961   -1.7661   -1.0906 H       1  BUTANE      0.0263
        12 H           4.9760    0.2731   -1.2896 H       1  BUTANE      0.0230
        13 H           4.9499   -1.3426   -2.0168 H       1  BUTANE      0.0230
        14 H           4.9778   -1.1679   -0.2537 H       1  BUTANE      0.0230
    @<TRIPOS>BOND
         1     1     2    1
         2     2     3    1
         3     3     4    1
         4     1     5    1
         5     1     6    1
         6     1     7    1
         7     2     8    1
         8     2     9    1
         9     3    10    1
        10     3    11    1
        11     4    12    1
        12     4    13    1
        13     4    14    1

No matter how you generate a mol2 file, if it corresponds to a molecule that can react, you must edit the mol2 file to give its reactive atoms **unique names**.  These atoms are referred to in the configuration file by name.

Configuration files
^^^^^^^^^^^^^^^^^^^

An HTPolyNet configuration file is a list of keyword:value pairs in YAML format.  The Library subpackage as a few example configuration files in the cfg directory.

The items in a configuration file break down into two major classes:

1. Items that specify run-time *parameters*.
2. Items that specify the system *chemistry*.

HTPolyNet parameters
''''''''''''''''''''

Below is a table of parameter keywords and descriptions for parameters that govern the overall execution of HTPolyNet.

General parameters:

===============================    ==============  =====================
Parameter                          Type            Description (default)
===============================    ==============  =====================
``Title``                          str             Descriptive title
``gmx_options``                    quoted string   options to pass to all ``gmx`` calls ('-quiet -nobackup')
``initial_density``                float           initial system density in kg/m3 (300.0)
``CURE_initial_search_radius``     float           initial capture radius in nm (0.5)
``CURE_radial_increment``          float           capture radius increment (0.25)
``desired_conversion``             float           desired fraction of possible crosslink bonds to form
``max_CURE_iterations``            int             maximum number of CURE iterations to run prior to reaching desired conversion
``late_treshhold``                 float           conversion above which reactions are all treated with probability 1.0
``charge_method``                  string          "gas" for Gasteiger; "bcc" for bcc; passed directly to antechamber
===============================    ==============  =====================

Other control parameters govern detailed aspects of the CURE algorithm.  These involve MD simulations performed immediately prior to and immediately after new bond addition to the topology, in order to relax those bonds.

Prior to introducing new bonds, one has the option of *"dragging"* atoms destined to be bonded to each other closer together in a series of dragging simulations.  The series is composed of stages, each of which involves three ``gmx mdrun`` calls: (1) a minimization; (2) an NVT relaxation; and (3) an NPT relaxation.  Soon-to-be-bonded atoms are connected by fictitious (type-6) harmonic bonds with equilibrium distances set at the current separation distances and relatively weak spring constants.  With each successive stage, the bond lengths are reduced and the spring constants increased until the desired separation distance and spring constant are achieved.  Dragging is optional.

Dragging parameters:

===============================    ==============  =====================
Parameter                          Type            Description (default)
===============================    ==============  =====================
``max_drag_stages``                int             number of drag stages to perform
``drag_limit``                     float           minimum distance each separation should achieve (nm); 0.0 turns off dragging (0.0)
``drag_nvt_steps``                 int             number of MD steps for NVT relaxation during dragging (-2, signals ``gmx mdrun`` to use the value in the mdp file)
``drag_npt_steps``                 int             number of MD steps for NPT relaxation during dragging (-2, signals ``gmx mdrun`` to use the value in the mdp file)
===============================    ==============  =====================

*After* new bonds are formed and all other bonded interactions, atom types, and charges are mapped from each bond's appropriate template, a series of *bond relaxation* MD simulations are performed.  These are in all ways similar to the optional *dragging* simulations except for the fact that here, the actual chemical bond parameters are progressively brought to their correct values as specified in the GAFF.  Bond relaxation is *required* because most new bonds are much longer than they should be at equilibrium.

Bond relaxation parameters:

===============================    ==============  =====================
Parameter                          Type            Description (default)
===============================    ==============  =====================
``max_bond_relaxation_stages``     int             number of bond relaxation stages to perform
``relax_nvt_steps``                int             number of MD steps for NVT relaxation 
``relax_npt_steps``                int             number of MD steps for NPT relaxation 
===============================    ==============  =====================

Chemistry parameters
''''''''''''''''''''

The system chemistries and initial composition are specified by a set of inter-referential YAML entries.

Top-level chemistry parameters:

================================= ====          ===========
Parameter                         Type          Description
================================= ====          ===========
``initial_composition``           dict          keys: monomer names, values: numbers of molecules in system
``reactions``                     list          reaction dicts, one per reaction
``use_symmetry_equivalent_atoms`` list          monomers for which symmetry-equivalent atoms are used
================================= ====          ===========

Reaction dicts:

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

Atom dicts:

=================== ====  ===========
Keyword             Type  Description
=================== ====  ===========
``reactant``        key   key to reactant in reactant dict to which this atom max_bond_relaxation_stages
``resid``           int   residue index in reactant molecular sequence to which this atom belongs (begins at 1)
``atom``            str   name of atom within that residue
``z``               int   number of available crosslink bonds for this atom
=================== ====  ===========

Bond dicts:

============= ======= ===========
Keyword       Type    Description
============= ======= ===========
``atoms``     list    list with the two atom keys the comprise the bond
``order``     float   bond order (currently not used; we let antechamber decide)
============= ======= ===========

