Usage
-----

Installation of the ``HTPolyNet`` package gives access to the ``htpolynet`` command:

.. code-block:: console

    $ htpolynet --help
    usage: htpolynet [-h] {run,parameterize,info,plots} ...

    HTPolyNet
    https://abramsgroup.github.io/HTPolyNet/

    Ming Huang
    mh3429@dragons.drexel.edu

    Cameron F. Abrams
    cfa22@drexel.edu

    Supported in part by Grants W911NF-17-2-0227 
    and W911NF-12-R-0011 from the US Army Research Lab

    positional arguments:
      {run,parameterize,info,plots}
        run                 build a system using instructions in the config file and any required molecular structure inputs
        parameterize        parameterize monomers and oligomer templates using instructinos in the config file
        info                print some information to the console
        plots               generate some plots that summarize aspects of the current completed build

    options:
      -h, --help            show this help message and exit

Subcommands
-----------

``htpolynet run``
^^^^^^^^^^^^^^^^^

This subcommand is for building crosslinked systems:

.. code-block:: console

  $ htpolynet run --help
  usage: htpolynet run [-h] [-lib LIB] [-proj PROJ] [-diag DIAG] [-restart] [--force-parameterization] [--force-checkin] [--loglevel LOGLEVEL] config

  positional arguments:
    config                input configuration file in YAML format

  options:
    -h, --help            show this help message and exit
    -lib LIB              local user library of molecular structures and parameterizations
    -proj PROJ            project directory; "next" (default) generates next directory 
                          Anything other than "next": 
                              if it exists, "-restart" must be included as a parameter; 
                              if not, it is created as a new project
    -diag DIAG            diagnostic log file
    -restart              restart in latest proj dir
    --force-parameterization
                          force GAFF parameterization of any input mol2 structures
    --force-checkin       force check-in of any generated parameter files to the system library
    --loglevel LOGLEVEL   Log level for messages written to diagnostic log (debug|info)

The arguments are explained in more detail below.

* ``config`` refers to the name of the :ref:`configuration input file <configuration_files>`.
* ``-lib`` names a directory that is treated as a local library of molecular structures.  By default, this is assumed to be ``./lib/``.  At the beginning of a new run, it should have one subdirectory ``molecules``.  Under ``molecules`` should be the two directories ``inputs`` and ``parameterized``.  ``HTPolyNet`` will look for input ``mol2`` or ``pdb`` files in ``molecules/inputs``, and put the results of parameterized molecules (i.e., Gromacs format ``gro``, ``itp``, and ``top`` files) in ``molecules/parameterized``.
* ``-proj`` names the project directory.  If provided with the value ``next`` (the default), ``HTPolyNet`` will create the next autonamed project directory.  These are always named as ``proj-``n, where n is replaced by an integer beginning with 0; if no project directory exists and an explicit one is not specified by ``-proj``, ``HTPolyNet`` creates the first one, ``proj-0``.  A project directory will automatically be given the following subdirectories:

  * ``molecules/parameterized`` -- all molecular parameterization results appear here
  * ``systems`` -- system initialization, CURE iterations, and post-cure all get their own subdirectories here.
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

The intention is that later invocations of ``htpolynet run`` can use these parameterizations without having to reperform them.  It is not required to parameterize outside of ``htpolynet run`` but sometimes it is useful to do so.

``htpolynet info``
^^^^^^^^^^^^^^^^^^

This simply outputs some information about ``HTPolyNet``.

``htpolynet plots``
^^^^^^^^^^^^^^^^^^^

If invoked inside of a directory containing one or more project directories, this instructs ``HTPolyNet`` to generate some plots.  This feature is currently under development.

