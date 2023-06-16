.. _configuration_postsim:

Configuration Files for ``htpolynet postsim``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``postsim`` is a subcommand that provides for convenient post-build MD simulations for annealing, equilibration, and making measurements of properties, like Young's modulus and glass transition temperature.  Its main advantage is that it understands the project file system structure and uses many sensible default values, making it easy to run post-build MD simulations.

A ``postsim`` configuration file is written in YAML and structured as a list of dictionaries, each entry of which is a *stage* which together constitute a *sequence* of stages.  Each stage involves an MD simulation for which you can specify the input coordinates and ``mdp`` parameters.  Because ``postsim`` executes a sequence of stages, the output of one stage can be the input for a later stage.

``HTPolyNet`` recoginizes the following types of post-build MD simulations

==============================  ===============  =================
Simulation Type                 keyword          Description
==============================  ===============  =================
Equilibration                   ``equilibrate``  Plain MD in NVT or NPT
Annealing                       ``anneal``       Temperature annealing using a low T and a high T and a schedule for cycling between them
Temperature-ladder              ``ladder``       Controlled change in T from one temperature to another over some number of stages of some duration
Uniaxial deformation            ``deform``       Uniaxial strain at a constant rate in one direction for a specified duration
==============================  ===============  =================

* Equilibration directives

  An ``equilibrate`` directive allows specification of several parameters:

  ================== ===========
  Keyword            Description
  ================== ===========
  ``subdir``         Name of directory in which simulation is conducted
                     relative to the project directory
                     default: ``postsim/equilibrate``
  ``input_top``      Input topology file relative to the project directory;
                     default: ``systems/final-results/final.top``
  ``input_gro``      Input coordinate file relative to the project directory;
                     default: ``systems/final-results/final.gro``
  ``input_grx``      Input extended atom attribute file relative to the
                     project directory; default: ``systems/final-results/final.grx``
  ``output_deffnm``  default gromacs output file basename;
                     default: ``equilibrate``
  ``traces``         List of trace types to generate output files for
  ``scatter``        Tuple specifying a scatter plot to generate based on the traces
  ``gromacs``        Gromacs directives; see :ref:`here <gromacs_yaml_directive>`         
  ``T``              Temperature in K
                     default: 300
  ``P``              Pressure in bar
                     default: 1
  ``ps``             Duration of simulation in picoseconds;
                     default: 1000
  ================== ===========

  A minimal post-build equilibration file for 100 ns at 310 K and 5 bar would contain:

  .. code-block:: YAML

    - equilibrate:
        ps: 100000
        T: 310
        P: 5

  Here, it uses the topology and coordinate files at the end of the build, and runs for 100 ns, with the results in ``<proj-dir>/postsim/equilibrate``.  Notice the ``-`` character at the left-most indent level; this indicates that the file is a list of items; the indentation of ``ps:`` indicates it that the identifier ``equilibrate`` is a keyword whose value is a dictionary.

  * Traces 
      
    ``postsim`` provides a shortcut approach to using ``gmx energy`` to extract energy-like quantities from Gromacs energy file (``edr``).  For example, if we would like to extract traces of temperature and density:

    .. code-block:: YAML

      - equilibrate:
          ps: 1000
          traces: 
            - Temperature
            - Density

    Note we have to list the names of the energy-like quantities exactly how ``gmx energy`` reports them.

  * Scatter plots 
    
    ``postsim`` also provides a shortcut interface to ``matplotlib`` to generate simple trace plots of any quantity extracted using ``traces``.  The single argument of the ``scatter`` keyword is a tuple in which the first element is the name of the x-data (typically ``time(ps)``), the second is a list of names of energy-like quantities to trace (e.g., ``[Density, Temperature]``), and the third element is the name of an output image file to contain the plot.

    We can add a ``scatter`` directive thusly

    .. code-block:: YAML

      - equilibrate:
          ps: 1000
          T: 310
          P: 5
          traces:
            - Temperature
            - Density
          scatter:
            - time(ps)
            - - Density 
              - Temperature
            - my_plot.png

  * Running ``postsim``

    Suppose we put the YAML shown above into the file ``eq.yml``.  To use this input configuration to launch a new post-build equilibration in the project directory ``proj-1``, one would issue the command

    .. code-block:: bash

      htpolynet postsim -cfg eq.yml -proj proj-1 -ocfg <build.cfg>

    The ``-ocfg`` key identifies the YAML config file used for the actual build; ``postsim`` extracts some information about the monomer library and other inputs from here; it also uses the Gromacs configuration of the original build config unless overridded in the postsim config.

* Annealing

  Annealing simulations are useful to help in densifying amorphous systems; the idea is that cycles of heating and cooling help the atoms fall deeper into local potential energy wells.  It is usually a good idea to anneal a newly built system before measuring any properties.

  An ``anneal`` directive allows specification of several parameters:

  ================== ===========
  Keyword            Description
  ================== ===========
  ``subdir``         Name of directory in which simulation is conducted
                     relative to the project directory
                     default: ``postsim/anneal``
  ``input_top``      Input topology file relative to the project directory;
                     default: ``systems/final-results/final.top``
  ``input_gro``      Input coordinate file relative to the project directory;
                     default: ``systems/final-results/final.gro``
  ``input_grx``      Input extended atom attribute file relative to the
                     project directory; default: ``systems/final-results/final.grx``
  ``output_deffnm``  default gromacs output file basename;
                     default: ``equilibrate``
  ``traces``         List of trace types to generate output files for
  ``scatter``        Tuple specifying a scatter plot to generate based on the traces
  ``gromacs``        Gromacs directives; see :ref:`here <gromacs_yaml_directive>`
  ``P``              Pressure (bar)
  ``ncycles``        number of annealing cycles
  ``T0``             "low" temperature
  ``T1``             "high" temperature
  ``T0_to_T1_ps``    duration of T-change simulation for T0 to T1 in ps
  ``T1_ps``          duration of T-hold at T1 in ps
  ``T1_to_T0_ps``    duration of T-change simulaiton for T1 to T0 in ps
  ``T0_ps``          duration of T-hold at T0 in ps (completing a cycle)
  ================== ===========

  An example ``anneal`` postsim directive might look like:

  .. code-block:: YAML

    - anneal:
        input_top: 'systems/final-results/final.top'
        input_gro: 'systems/final-results/final.gro'
        P: 1
        T0: 300
        T1: 600
        ncycles: 3
        T0_to_T1_ps: 10000
        T1_ps: 10000
        T1_to_T0_ps: 10000
        T0_ps: 10000


* Temperature ladders

  Temperature-ladder simulations are useful in the response of system density to temperature over a wide temperature domain.  This can be useful for estimating the coefficient of thermal expansion and the glass transition temperature.  

  A ``ladder`` directive allows specification of several parameters:

  ================== ===========
  Keyword            Description
  ================== ===========
  ``subdir``         Name of directory in which simulation is conducted
                     relative to the project directory
                     default: ``postsim/anneal``
  ``input_top``      Input topology file relative to the project directory;
                     default: ``systems/final-results/final.top``
  ``input_gro``      Input coordinate file relative to the project directory;
                     default: ``systems/final-results/final.gro``
  ``input_grx``      Input extended atom attribute file relative to the
                     project directory; default: ``systems/final-results/final.grx``
  ``output_deffnm``  default gromacs output file basename;
                     default: ``equilibrate``
  ``traces``         List of trace types to generate output files for
  ``scatter``        Tuple specifying a scatter plot to generate based on the traces
  ``gromacs``        Gromacs directives; see :ref:`here <gromacs_yaml_directive>`
  ``P``              Pressure (bar)
  ``Tlo``            "low" temperature (K)
  ``Thi``            "high" temperature (K)
  ``deltaT``         change in temperature between subsequent simulations on the ladder
  ``ps_per_rise``    duration (ps) of T-change simulation beginning a new step on the ladder
  ``ps_per_run``     duration (ps) of T-hold simulation at each step
  ``warmup_ps``      duration (ps) of a warmup simulation before beginning T-changes
  ================== ===========

  An example ``ladder`` postsim directive might look like:

  .. code-block:: YAML

    - ladder:
        input_top: 'systems/final-results/final.top'
        input_gro: 'postsim/equilibrate/equilibrate.gro'
        subdir: 'postsim/ladder-heat'
        Tlo: 300
        Thi: 600
        deltaT: 5
        ps_per_rise: 100
        ps_per_run: 900
        warmup_ps: 1000 

  One can go either up or down a ladder.  If ``deltaT`` is specified as a postive number, the initial temperature is ``Tlo``.  However, if ``deltaT`` is specified as a negative number, the intial temperature is ``Thi`` instead of ``Tlo``.  One way I like to use ladders is as one big cycle, measuring thermal response in both directions:

  .. code-block:: YAML

    - ladder:
        input_top: 'systems/final-results/final.top'
        input_gro: 'postsim/equilibrate/equilibrate.gro'
        subdir: 'postsim/ladder-heat'
        Tlo: 300
        Thi: 600
        deltaT: 5
        ps_per_rise: 100
        ps_per_run: 900
        warmup_ps: 1000 
    - ladder:
        input_top: 'systems/final-results/final.top'
        input_gro: 'postsim/ladder-heat/ladder.gro'
        subdir: 'postsim/ladder-cool'
        Tlo: 300
        Thi: 600
        deltaT: -5
        ps_per_rise: 100
        ps_per_run: 900
        warmup_ps: 1000 

  A measure of how well equilibrated your system is and how carefully you conducted the ladder simulations is agreement between the thermal responses in the T-up and T-down ladders.

  By default, a ``ladder`` generates a trace of both temperature and density from the associated ``edr`` files.  The ``plots`` subcommand can be used to read these files for a set of equivalent project directories and compute a glass transition temperature.

* Uniaxial deformation

  Constant strain-rate uniaxial deformation MD simulations can be performed using ``postsim`` as well.  Such simulations are used to generate stress-strain curves from which moduli are extracted.  

  A ``deform`` directive allows specification of several parameters:

  ================== ===========
  Keyword            Description
  ================== ===========
  ``subdir``         Name of directory in which simulation is conducted
                     relative to the project directory
                     default: ``postsim/anneal``
  ``input_top``      Input topology file relative to the project directory;
                     default: ``systems/final-results/final.top``
  ``input_gro``      Input coordinate file relative to the project directory;
                     default: ``systems/final-results/final.gro``
  ``input_grx``      Input extended atom attribute file relative to the
                     project directory; default: ``systems/final-results/final.grx``
  ``output_deffnm``  default gromacs output file basename;
                     default: ``equilibrate``
  ``traces``         List of trace types to generate output files for
  ``scatter``        Tuple specifying a scatter plot to generate based on the traces
  ``gromacs``        Gromacs directives; see :ref:`here <gromacs_yaml_directive>`
  ``direction``      Direction along with strain is applied; x, y, or z
  ``edot``           Dimensionless strain rate in ps\ :sup:`-1`\; 0.0001 = 10\ :sup:`8`\  s\ :sup:`-1`\
  ``ps``             Duration of the simulation in ps; product of ``ps`` and ``edot`` gives the strain at the end of the simulation.
  ================== ===========

  An example ``deform`` postsim directive might look like:

  .. code-block:: YAML

    - deform:
        input_top: 'systems/final-results/final.top'
        input_gro: 'postsim/equilibrate/equilibrate.gro'
        subdir: 'postsim/deform-x'
        direction: x
        edot: 0.0001
        ps: 1000  # final strain 10%

  For isotropic (bulk) systems, the x, y, and z directions are equivalent, and one should measure material response along all three directions.  One can use a single input configuration file to stipulate all three such simulations, each starting from a common initial condition:

  .. code-block:: YAML

    - deform:
        input_top: 'systems/final-results/final.top'
        input_gro: 'postsim/equilibrate/equilibrate.gro'
        subdir: 'postsim/deform-x'
        direction: x
        edot: 0.0001
        ps: 1000  # final strain 10%
    - deform:
        input_top: 'systems/final-results/final.top'
        input_gro: 'postsim/equilibrate/equilibrate.gro'
        subdir: 'postsim/deform-y'
        direction: y
        edot: 0.0001
        ps: 1000  # final strain 10%
    - deform:
        input_top: 'systems/final-results/final.top'
        input_gro: 'postsim/equilibrate/equilibrate.gro'
        subdir: 'postsim/deform-z'
        direction: z
        edot: 0.0001
        ps: 1000  # final strain 10%

  By default, a ``deform`` block generates a traces of normal stress tensor components along the pulling direction and box size along that direction.  This allows the ``plots`` subcommand to estimate Young's modulus.
