Usage
-----

Installation of the ``HTPolyNet`` package gives access to the ``htpolynet`` command.  Invoking ``htpolynet --help`` shows basic usage of the command, which always involves providing a *subcommand*:

.. code-block:: console

    $ htpolynet --help
    usage: htpolynet [-h] {run,parameterize,info,plots,fetch-example,input-check} ...

    HTPolyNet 0.0.1
    https://abramsgroup.github.io/HTPolyNet/

    Ming Huang
    mh3429@dragons.drexel.edu

    Cameron F. Abrams
    cfa22@drexel.edu

    Supported in part by Grants W911NF-17-2-0227 
    and W911NF-12-R-0011 from the US Army Research Lab

    positional arguments:
      {run,parameterize,info,plots,fetch-example,input-check}
        run                 build a system using instructions in the config file and any 
                            required molecular structure inputs
        parameterize        parameterize monomers and oligomer templates using instructinos 
                            in the config file
        info                print some information to the console
        plots               generate some plots that summarize aspects of the current 
                            completed build
        fetch-example       fetch and unpack example(s) from the HTPolyNet.Library:
                           "0-liquid-styrene", "1-polystyrene", "2-polymethylstyrene",
                            "3-bisgma-styrene-thermoset", "4-pacm-dgeba-epoxy-thermoset", 
                            "5-dfda-fde-epoxy-thermoset", "6-htpb-ipdi"
        input-check         reports number of atoms that would be in initial system based 
                            on config

    options:
      -h, --help            show this help message and exit


``htpolynet run``
^^^^^^^^^^^^^^^^^

This subcommand is used to build a polymerized system given directives in the configuration file and any required input monomer structures.


The arguments for ``htpolynet run`` are explained in more detail below.

* ``config`` refers to the name of the :ref:`configuration input file <configuration_files>`.
* ``-lib`` names a directory that is treated as a local library of molecular structures.  By default, this is assumed to be ``./lib/`` (that is, the directory you are in when you issue ``htpolynet run`` is expected by default to have a ``lib/`` directory; if not, the Library subpackage of the ``HTPolyNet`` package will be queried for any data).  At the beginning of a new run, ``lib/`` should have one subdirectory called ``molecules``.  Under ``molecules`` should be the two directories ``inputs`` and ``parameterized``.  ``HTPolyNet`` will look for input ``mol2`` or ``pdb`` files in ``lib/molecules/inputs``, and "check-in" the results of parameterized molecules (i.e., Gromacs format ``gro``, ``itp``, and ``top`` files) in ``lib/molecules/parameterized``.
* ``-proj`` names the project directory.  If provided with the value ``next`` (the default), ``HTPolyNet`` will create the next autonamed project directory.  These are always named as ``proj-``n, where n is replaced by an integer beginning with 0; if no project directory exists and an explicit one is not specified by ``-proj``, ``HTPolyNet`` creates the first one, ``proj-0``.  A project directory will automatically be given the following subdirectories:

  * ``molecules/parameterized`` -- all molecular parameterization results appear here (in addition to being checked in to the library)
  * ``systems`` -- system initializations, equilibrations, CURE iterations, and postcure equilibrations all get their own subdirectories here.
  * ``plots`` -- various plots generated on the fly.

  These will be explained more fully in the tutorials.

* ``-diag`` names the diagnostic output file, and ``--loglevel`` specifies the logging level it uses.  The default level is ``debug`` (i.e., the most informative).
* ``-restart`` indicates that this is a restart (experimental!).
* ``--force-parameterization`` signals that ``HTPolyNet`` should perform all molecular parameterizations from scratch even if parameterizations exist in the library.
* ``--force-checkin`` signals that any parameterizations ``HTPolyNet`` does should have their results "checked-in" to the library, even if previous parameterizations are there already.

``htpolynet parameterize``
^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the command for only performing molecular parameterizations and checking the results into a library:

.. code-block:: console

  $ htpolynet parameterize --help
  usage: htpolynet parameterize [-h] [-lib LIB] [-diag DIAG] [-restart] [--force-parameterization] [--force-checkin] [--loglevel LOGLEVEL] config

  positional arguments:
    config                input configuration file in YAML format

  options:
    -h, --help            show this help message and exit
    -lib LIB              local user library of molecular structures and parameterizations
    -diag DIAG            diagnostic log file
    -restart              restart in latest proj dir
    --force-parameterization
                          force GAFF parameterization of any input mol2 structures
    --force-checkin       force check-in of any generated parameter files to the system library
    --loglevel LOGLEVEL   Log level for messages written to diagnostic log (debug|info)

The command-line options of ``htpolynet parameterize`` have all the same meanings as they do for ``htpolynet run``.  The only difference is that ``htpolynet parameterize`` **only** performs the parameterization of all monomers and oligomer templates.  The intention is that later invocations of ``htpolynet run`` can use these parameterizations without having to reperform them.  Of course, since a first invocation of ``htpolynet run`` also peforms parameterizations and saves their results, it is strictly **never** necessary to use ``htpolynet parameterize``.  However, if your parameterizations have issues, it is cleaner to use ``htpolynet parameterize`` to try to fix them.

``htpolynet info``
^^^^^^^^^^^^^^^^^^

This simply outputs some information about ``HTPolyNet``.

.. code-block:: console

  $ htpolynet info
  This is some information on your installed version of HTPolyNet
  System library is /home/cfa/Git/HTPolyNet/Library
  Ambertools:
   antechamber (ver. 22.0) at antechamber                                       
         tleap (ver. 22.0) at tleap                                             
      parmchk2 (ver. 22.0) at parmchk2 

``htpolynet info`` only reports the absolute pathname of the ``Library`` subpackage for your reference, and the fully resolved command names for the three required Ambertools executables ``antechamber``, ``tleap`` and ``parmchk2``, along with their versions.  If they are already in your path, the results appear as above. 

``htpolynet plots``
^^^^^^^^^^^^^^^^^^^

If invoked inside of a directory containing one or more project directories, this instructs ``HTPolyNet`` to generate some plots.

.. code-block:: console

  $ htpolynet plots --help
  usage: htpolynet plots [-h] [-logs LOGS [LOGS ...]] [-proj PROJ] [-o O] [--plotfile PLOTFILE]

  options:
    -h, --help            show this help message and exit
    -logs LOGS [LOGS ...]
                          names of diagnostic log files (1 or more)
    -proj PROJ            name of project directory
    -o O                  name of global trace output data file
    --plotfile PLOTFILE   name of plot file to generate

We explain detailed usage of ``htpolynet plots`` in the tutorials.  Briefly, if a project directory is name via the ``-proj`` option, ``HTPolyNet`` will generate a set of plots tracing the system temperature, density, and number of polymerization bonds vs simulation time.  If one or more diagnostic log files is named in the ``-logs`` option, ``HTPolyNet`` will generate a pair of plots of conversion vs. wall-clock time and iteration vs wall-clock time, with each diagnostic log getting its own curve on each plot.

``htpolynet fetch-example``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This will fetch one or more examples from the ``example_depot`` of the system ``Library``.

.. code-block:: console

  $ htpolynet fetch-example --help
  usage: htpolynet fetch-example [-h] [-n {0,1,2,3,4,5,6,all}] [-k]

  options:
    -h, --help            show this help message and exit
    -n {0,1,2,3,4,5,6,all}
                          number of example tarball to unpack from 0-liquid-styrene, 1-polystyrene,
                          2-polymethylstyrene, 3-bisgma-styrene-thermoset, 4-pacm-dgeba-epoxy-thermoset, 
                          5-dfda-fde-epoxy-thermoset,
                          6-htpb-ipdi
    -k                    keep tarballs

Fetching will copy the tarball for the requested system to the current directory and then untar it and remove it, leaving behind the directory.  For example to fetch the PACM-DGEBA epoxy thermoset example:

.. code-block:: console

  $ htpolynet fetch-example -n 4
  $ ls
  4-pacm-dgeba-epoxy-thermoset/
  $ cd 4-pacm-dgeba-epoxy-thermoset
  $ ls
  DGE-PAC-hi.yaml  DGE-PAC-lo.yaml  lib/  run.sh

This folder (like all example folders) comes with two configuration files that differ only the the requested degree of cure.  "hi" refers to 95\% cure, and "lo" to 50\%.  Also provided is the ``./lib/molecules`` folders with the ``./lib/molecules/inputs`` and ``./lib/molecules/parameterized`` empty subfolders.  Finally, the bash script ``run.sh`` can just be invoked to build the input monomers and run the two builds in series.  This will be described in much more detail in the tutorials.

``htpolynet fetch-example -n all`` just grabs all seven examples.

``htpolynet input-check``
^^^^^^^^^^^^^^^^^^^^^^^^^

The purpose of this subcommand is to report the size of the initial system that *would* be created by the provided configuration file and monomer input structures.

.. code-block:: console

  $ htpolynet input-check DGE-PAC-hi.yaml
  Molecule DGE: 53 atoms, 200 molecules
  Molecule PAC: 41 atoms, 100 molecules
  DGE-PAC-hi.yaml: 14700 atoms in initial system
