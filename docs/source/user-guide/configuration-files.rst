.. _configuration_files:

Configuration Files
-------------------

The configuration file is how the user tells ``HTPolyNet`` what it needs in order to generate a polymerized system beginning with structures of the individual monomers and a description of the polymerization chemistry.  ``HTPolyNet`` expects configuration files to be in ``YAML`` format.  In essence, the configuration is a set of dictonaries of key:value pairs, with each dictionary corresponding to certain phases of ``HTPolyNet`` execution.

For example, a simple configuration file that describes building a system of polystyrene from a liquid of styrene monomers might look like::

    Title: polystyrene
    gromacs: {
        gmx: 'gmx',
        gmx_options: '-quiet -nobackup',
        mdrun: 'gmx mdrun'
    }
    ambertools: {
        charge_method: gas
    }
    constituents: {
        STY: {count: 100}
    }
    densification: {
        initial_density: 300.0,  # kg/m3
        equilibration: [
            { ensemble: min },
            { ensemble: nvt, temperature: 300, ps: 10 },
            { ensemble: npt, temperature: 300, pressure: 10, ps: 200 }
        ]
    }
    precure: {
        preequilibration: {
            ensemble: npt,
            temperature: 300,        # K
            pressure: 1,             # bar
            ps: 200
        },
        anneal: {
            ncycles: 2,
            initial_temperature: 300,
            cycle_segments: [
                { T: 300, ps: 0 },
                { T: 600, ps: 20 },
                { T: 600, ps: 20 },
                { T: 300, ps: 20 },
                { T: 300, ps: 20 }
            ]
        },
        postequilibration: {
            ensemble: npt,
            temperature: 300,        # K
            pressure: 1,             # bar
            ps: 100
        }
    }
    CURE: {
        controls: {
            initial_search_radius: 0.5, # nm
            radial_increment: 0.25,     # nm
            max_iterations: 150, 
            desired_conversion: 0.95,
            late_threshhold: 0.85
        },
        drag: {
            trigger_distance: 0.6,   # nm
            increment: 0.08,         # nm
            limit: 0.3,              # nm
            equilibration: [
                { ensemble: min },
                { ensemble: nvt, temperature: 600, nsteps: 1000 },
                { ensemble: npt, temperature: 600, pressure: 1, nsteps: 2000 }
            ]
        },
        relax: {
            increment: 0.08,         # nm
            equilibration: [
                { ensemble: min },
                { ensemble: nvt, temperature: 600, nsteps: 1000 },
                { ensemble: npt, temperature: 600, pressure: 1, nsteps: 2000 }
            ]
        },
        equilibrate: {
            ensemble: npt,
            temperature: 300,       # K
            pressure: 1,            # bar
            ps: 100
        },
        gromacs: {
            rdefault: 0.9 # nm
        }
    }
    postcure: {
        anneal: {
            ncycles: 2,
            initial_temperature: 300,
            cycle_segments: [
                { T: 300, ps: 0 },
                { T: 600, ps: 20 },
                { T: 600, ps: 20 },
                { T: 300, ps: 20 },
                { T: 300, ps: 20 }
            ]
        },
        postequilibration: {
            ensemble: npt,
            temperature: 300,       # K
            pressure: 1,            # bar
            ps:  100
        }
    }
    reactions:
    - {
        name:        'sty1_1',
        stage:       cure,
        reactants:   {1: STY, 2: STY},
        product:     STY~C1-C2~STY,
        probability: 1.0,
        atoms: {
            A: {reactant: 1, resid: 1, atom: C1, z: 1},
            B: {reactant: 2, resid: 1, atom: C2, z: 1}
        },
        bonds: [
            {atoms: [A, B], order: 1}
        ]
      }
    - {
        name:         'styCC',
        stage:        cap,
        reactants:    {1: STY},
        product:      STYCC,
        probability:  1.0,
        atoms: {
            A: {reactant: 1, resid: 1, atom: C1, z: 1},
            B: {reactant: 1, resid: 1, atom: C2, z: 1}
        },
        bonds: [
            {atoms: [A, B], order: 2}
        ]
      }

This example file contains nine distinct **directives**.  

* ``gromacs``:  This directive specifies parameters ``HTPolyNet`` uses when invoking the Gromacs executable.
  
    =====================================    ==============  =====================
    ``gromacs`` subdirective                 Type            Description (default)
    =====================================    ==============  =====================
    ``gmx``                                  str             ``gmx`` command (default ``gmx``)
    ``gmx_options``                          quoted string   options to pass to all ``gmx`` calls (default ``-quiet -nobackup``)
    ``mdrun``                                quoted string   ``mdrun`` command (default ``gmx (options) mdrun``)
    ``mdrun_single_molecule``                quoted string   version of ``mdrun`` to use for any single-molecule Gromacs runs
    ``mdrun_options``                        dict            command-line arguments to pass to ``mdrun`` (none)
    =====================================    ==============  =====================

    If you are running on a supercomputer with a native installation of Gromacs, it is likely you should point the parameter ``gmx`` to the fully resolved pathname of ``gmx_mpi`` (or load the appropriate module), and use the ``mdrun`` parameters to specify the ``mpirun`` or ``mpiexec`` syntax needed to launch ``gmx_mpi mdrun``.  The ``gromacs_single_molecule`` subdirective allows you to specify a particular form of ``mdrun`` appropriate for single-molecule simulations.  These are most often used as part of parameterization or conformer generation.  Typically, it's best to run these on a single processor without domain decomposition.

    The ``gromacs`` directive is optional; if none is specified the default values are used.

* ``ambertools``:  This directive specifies parameters ``HTPolyNet`` uses when working with the AmberTools suite.

    =====================================    ==============  =====================
    ``ambertools`` subdirective              Type            Description (default)
    =====================================    ==============  =====================
    ``charge_method``                        string          charge model used by ``antechamber`` (default ``gas``)
    =====================================    ==============  =====================

    For now, you can choose any charging method compatible with ``antechamber``.  The ``antechamber`` directive is optional.

* ``constituents``
  
    This **required** directive is a set of one or more "key":"record" pairs in which each key is the name of a molecule (here, "STY") and the record is a dictionary of keyword:value pairs.  The allowable keywords in a ``constituent`` record are as follows.

    =====================================    =================  =====================
    ``constituents`` record keyword          Type               Description (default)
    =====================================    =================  =====================
    ``count``                                int                (required) number of these molecules in the system
    ``stereocenters``                        list               (optional) list of names of chiral carbon atoms if any
    ``symmetry_equivalent_atoms``            list               (optional) list of sets of symmetry equivalent atom names, if any
    ``conformers``                           conformers record  (optional) parameters specifying if and how initial conformers are generated
    =====================================    =================  =====================

    In the example above, we are requesting a system of 100 styrene molecules.  The key ``STY`` signals to ``HTPolyNet`` that it should look for either ``STY.mol2`` or ``STY.pdb`` in ``./lib/molecules/inputs`` **or** it should look for ``STY.gro``, ``STY.itp``, ``STY.top``, and ``STY.grx`` in ``./lib/molecules/parameterized``.  The latter is the case if either ``htpolynet run`` or ``htpolynet parameterized`` has already been run with ``STY.mol2`` or ``STY.pdb``.  Multiple records in ``constituents`` should all have the "key":"record" syntax and be separated by commas.

    ``HTPolyNet`` allows you to use multiple conformers of flexible molecules to build the initial liquid system.  It can use either ``obabel``'s ``confomers`` capability or an MD simulation via ``gromacs`` to generate these.  The ``conformers`` record has two subdirectives:
    =====================================    ==========================  =====================
    ``conformers`` record keyword            Type                        Description (default)
    =====================================    ==========================  =====================
    ``count``                                int                         (required) number of unique conformers to generate (per stereoisomer)
    ``generator``                            conformer generator record  (optional) parameters specifying how conformers are generated
    =====================================    ==========================  =====================

    The ``conformers.generator`` record has several subdirectives:
    =====================================    ===========================  =====================
    ``generator`` record keyword             Type                         Description (default)
    =====================================    ===========================  =====================
    ``name``                                 str                          (required) ``obabel`` or ``gromacs``
    ``params``                               generator parameters record  (optional) parameters specifying the generator's operation (only relevant for ``gromacs``)
    =====================================    ===========================  =====================

    The ``conformers.generator.params`` record has several subdirectives:
    =====================================    ===========================  =====================
    ``params`` record keyword                Type                         Description (default)
    =====================================    ===========================  =====================
    ``ensemble``                             str                          ``nvt`` is the only option that makes sense
    ``temperature``                          float                        (optional) Temperature of the conformer-generating MD simulation
    ``ps``                                   float                        (optional) Duration of the conformer-generating MD simulation
    ``pad``                                  float                        (optional) Box-size padding for the vacuum MD simulation
    =====================================    ===========================  =====================

* ``densification``

    This directive instructs ``HTPolyNet`` how to run the initial densification of the fresh simulation system.  It has two subdirectives:

    =====================================    ==============  =====================
    ``densification`` subdirective           Type            Description (default)
    =====================================    ==============  =====================
    ``initial_density``                      float           density in kg/m^3 at which molecules are placed randomly into a box to make the initial coordinates (default 300.0)
    ``equilibration``                        list            list of **MD records** 
    =====================================    ==============  =====================

    The ``equilibration`` subdirective should contain one or more *MD records*. An MD record is a dictionary of keyword:value pairs:

    =====================================    ==============  =====================
    MD record keyword                        Type            Description
    =====================================    ==============  =====================
    ``ensemble``                             string          (required) min (minimization), npt, or nvt
    ``temperature``                          float           (required if ``ensemble`` is nvt or npt) Temperature in K assigned to ``ref_t`` in Gromacs ``mdp`` file
    ``pressure``                             float           (required if ``ensemble`` is npt) Pressure in bar assigned to ``ref_p`` in Gromacs ``mdp`` file
    ``nsteps``                               int             (optional; required if ``ps`` not provided) Duration of MD simulation in number of time steps
    ``ps``                                   float           (optional; required if ``nsteps`` not set) Duration of MD simulation in picoseconds
    ``repeat``                               int             (optional) number of times to repeat this simulation in series; default is 0 (i.e., run once)
    =====================================    ==============  =====================

    The ``repeat`` subdirective is especially useful for densifications that start at very low initial densities.  It is better to run several short NPT simulations than a single long one so that the box size shrinkage doesn't overwhelm Gromacs' domain decomposition algorithm.

* ``precure``
    
    The ``precure`` directive instructs ``HTPolyNet`` on running a series of MD simulations after densification but before the cure.  There are three allowable subdirectives for ``precure``: 

    =====================================    =================    =====================
    ``precure`` subdirective                 Type                 Description (default)
    =====================================    =================    =====================
    ``preequilibration``                     MD record            optional MD simulation
    ``anneal``                               **Anneal record**    Description of an annealing simulation after the optional ``preequilibration``
    ``postequilibration``                    MD record            optional MD simulation         
    =====================================    =================    =====================

    Both the ``preequilibration`` and ``postequilibration`` directives contain MD records described above.  The *Anneal record* has the following subdirectives:

    =====================================    =================    =====================
    Anneal record subdirective               Type                 Description (default)
    =====================================    =================    =====================
    ``ncycles``                              int                  number of annealing cycles
    ``initial_temperature``                  float                (optional) Initial temperaure in K, really only sets the ``gen-temp`` ``mdp`` parameter 
    ``cycle_segments``                       list                 list of **cycle records**
    =====================================    =================    =====================

    A **cycle record** corresponds to an "annealing-point" in the Gromacs ``mdp`` file.  

    =====================================    =================    =====================
    Cycle record subdirective                Type                 Description (default)
    =====================================    =================    =====================
    ``T``                                    float                Targe temperature in K 
    ``ps``                                   float                cycle duration; if prior ``T`` is different, simulation is *brought to* this ``T`` in this amount of time; if prior ``T`` is the same, simulation is *held at* this ``T`` for this amount of time.
    =====================================    =================    =====================

    Each cycle consists of one pass through the cycle segments.  In the example here, one cycle consists of Gromacs taking the system from 300 to 600 K in the first 20 ps, then holding at 600 for 20 pm, then reducing to 300 K over 20 ps and holding it there for 20 ps.

* ``CURE``
   
    This directive contains all instructions governing the :ref:`CURE algorithm <cure_section>`.  There are five possible subdirectives:

    =====================================    =================    =====================
    ``CURE`` subdirective                    Type                 Description (default)
    =====================================    =================    =====================
    ``controls``                             list                 Control parameter values
    ``drag``                                 list                 Dragging parameter values
    ``relax``                                list                 Bond relaxation parameter values
    ``equilibrate``                          MD record            CURE iteration equilibration parameters
    ``gromacs``                              list                 any ``mdp`` keyword:value pairs to include in all ``mdp`` files in the ``CURE`` sequence
    =====================================    =================    =====================

    * ``CURE.controls`` parameters

        =================================    =================   ======================
        ``CURE.controls`` parameter          Type                Description (default)
        =================================    =================   ======================
        ``initial_search_radius``            float               initial search radius in nm (default 0.5)
        ``radial_increment``                 float               increment by which search radius is increased if no bonds are found at current radius (default 0.25 nm)
        ``max_iterations``                   int                 absolute maximum number of allowed iterations (default 150), 
        ``desired_conversion``               float [0-1]         target conversion between 0 and 1.0 (default 0.95)
        ``late_threshhold``                  float [0-1]         conversion above which bond probabilities are ignored
        =================================    =================   ======================

.. _cure.drag:

    * ``CURE.drag`` parameters:  Dragging refers to a series of MD simulations (called "stages") in which harmonic restraints are applied to each pair of atoms assigned to form a bond, but **before** the bonds actually form.  Dragging is useful to reduce 1-4 distances that ultimately arise when bonds form.  Each stage in the series uses a specially modified topology file in which "new" bonds of type 6 are added, one for each pair of to-be-bonded atoms. Each of these bonds has a parameter ``kb``, the spring constant, and ``b0``, the equilibrium length.  The ``drag`` directive governs how those ``b0`` parameters are linearly decreased through the set of stages to slowly bring the atoms closer together.   The ``limit`` parameter is the target distance of dragging, and ``increment`` determines the number of stages it will take to get there.

        =================================    =================   ======================
        ``CURE.drag`` parameter              Type                Description (default)
        =================================    =================   ======================
        ``increment``                        float               minimum amount by which target ``drag`` distance is decreased in steps (default 0.08)
        ``limit``                            float               distance in nm to which all bonds are dragged (default 0.3)
        ``equilibration``                    MD record           describes the MD simulations used to equilibrate at each stage 
        =================================    =================   ======================

.. _cure.relax:

    * ``CURE.relax`` parameters:  Relaxation refers to a series of MD simulations (also called "stages") in which the ``kb`` and ``b0`` parameters of each new bond are "attenuated" from a weak (low ``kb``), long (large ``b0``) state to the state dictated by the force field.  The ``increment`` determines the number of stages are performed.

        =================================    =================   ======================
        ``CURE.relax`` parameter             Type                Description (default)
        =================================    =================   ======================
        ``increment``                        float               minimum amount by which ``b0`` parameters are decreased in steps (default 0.08)
        ``equilibration``                    MD record           describes the MD simulations used to equilibrate at each stage 
        =================================    =================   ======================

    * ``gromacs`` parameters:  These parameters govern modification to ``mdp`` files used in the dragging and relaxation MD simulations.  ``HTPolyNet`` adjusts the cutoff distances to conform to the longest unrelaxed bond in the system, and the ``rdefault`` parameter provides the floor below which it will not go any lower.

        =================================    =================   ======================
        ``CURE.gromacs`` parameter           Type                Description (default)
        =================================    =================   ======================
        ``rdefault``                         float               minimum cutoff radius (default 0.9)
        =================================    =================   ======================


* ``postcure`` 

    The ``postcure`` directive instructs ``HTPolyNet`` on running a series of MD simulations after cure.  Its form is identical to that of ``precure``, namely with optional ``preequilibration``, ``anneal``, and ``postequilibration`` subdirectives.

.. _reactions:

* ``reactions``

    The ``reactions`` directive contains a list of **reaction records**.  HTPolyNet expects one or more reaction templates to be defined in the configuration file.  A reaction is defined by the precise pairs of atoms that become new covalent bonds.  To precisely define each such pair, the reaction must also identify one or more reactant molecules.  Each reaction also names a single product molecule.  HTPolyNet will build oligomer templates using these reactions and then GAFF-parameterize them.  The parameterizations are used during CURE to re-type atoms and reset charges after each new bond is formed.

    ==============================  ==========  =================
    ``reaction`` record directives  Type        Description
    ==============================  ==========  =================
    ``name``                        str         descriptive name
    ``stage``                       str         one of ``cure``, ``cap``, ``build``, or ``param``
    ``probability``                 float       probability that bond will form in one iteration if identified (1.0)
    ``reactants``                   dict        keyword: reactant key, value: reactant molecule name
    ``product``                     str         name of product molecule
    ``atoms``                       dict        keyword: atom key, value: **atom record**
    ``bonds``                       list        list of **bond records**, one item per bond formed in reaction
    ==============================  ==========  =================

    The ``stage`` value signifies how ``HTPolyNet`` uses the reaction.  It will generate GAFF parameters and topologies for any product of a reaction with stage ``cure``, ``cap``, or ``param``.  ``cure`` reactions are those assigned to take place during CURE.  ``cap`` reactions are optional and take place once the CURE has finished; these can be used to revert the active form of any unreacted monomers back to their proper forms.  ``

    The ``atoms`` directive is a dictionary of atom records where the key is an atom "key", which is referenced in bond record.

    * Atom records uniquely identify atoms in reactants, assigning them a shorthand key that is used in subsequent bond records.
        
        ======================== ============== =================
        Atom record subdirective type           Description
        ======================== ============== =================
        ``reactant``             arb.           Reactant key that references the ``reactants`` directive of the reaction
        ``resid``                int            Residue number in the reactant containing this atom
        ``atom``                 str            Atom name (originates in monomer ``mol2`` or ``pdb`` file)
        ``z``                    int            Number of possible bonds atom can participate in
        ======================== ============== =================

    * Bond records specify the bond(s) that form during this reaction.

        ======================== ============== =================
        Bond record subdirective type           Description
        ======================== ============== =================
        ``atoms``                list           The two atom keys that define the atoms that form the bond
        ``order``                int            Order (1=single, 2=double) of resulting bond
        ======================== ============== =================

    In the example here, we define two unique reactions.  One is the C1-C2 bond that links two styrene monomers, and the other is the *intramolecular* C1-C2 double bond that "reverts" the active form of a monomer back to its "proper" form.  Since that reaction's ``stage`` is ``cap``, this signifies that it is formed only **after** CURE has finished.
