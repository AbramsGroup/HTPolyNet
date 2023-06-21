.. _configuration_analyze:

Configuration Files for ``htpolynet analyze``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``analyze`` is a subcommand that provides a shortcut interface to ``gmx`` subcommands for performing analyses of MD trajectories generated during ``run`` or ``postsim`` phases.

An ``analyze`` configuration file is structured as a list of dictionaries, each entry of which is a *stage* of the analysis.  Each stage consists of execution of a ``gmx`` analysis and parsing/reporting of its output.  

``htpolynet analyze`` is a *generic* way to interface with the ``gmx`` subcommands, and its main advantage is that it is aware of the HTPolyNet project directory structure.

A stage dictionary must have values for the ``command`` and ``subdir`` keywords; all others are optional.

================== ===========
Keyword            Description
================== ===========
``command``        Gromacs subcommand (required)
``subdir``         Name of directory in which analysis is conducted relative to the *project* directory (required)
``links``          List of files to make symlinks to in the analysis directory, with pathnames relative to the project directory
``gromacs``        Gromacs directives; see :ref:`here <gromacs_yaml_directive>`; default is the stipulation that the Gromacs command-line is ``gmx``
``options``        Keyword:value pairs for any options to the ``gmx`` subcommand.  Options with no values must be explicitly given an empty string value.
``outfile``        Name of text file to create and report any output of the command in.
``matchlines``     List of greppable strings; only lines matching one or more of these strings are output, if set
``console-input``  List of console inputs one would provide if running the ``gmx`` subcommand interactively, as strings
================== ===========

Simple examples
!!!!!!!!!!!!!!!

1.  ``gmx check``

  Suppose you wanted to run ``gmx check`` on ``proj-0/systems/densification/densified-npt.trr`` from the base directory.  At the command line, you could simply issue

  .. code-block:: bash

    $ gmx check -f proj-0/systems/densification/densified-npt.trr

  This is of course very simple.  If you want HTPolyNet to do the same thing, and put the output in a particular place, you could create the YAML input file ``my-analysis.yaml``:abbr:

  .. code-block:: yaml

    - command: check
      subdir: analyze/check
      links:
        - systems/densification/densified-npt.trr
      options:
        f: densified-npt.trr
      outfile: check.out
      matchlines: ['Last frame','Coords']

  This will generate ``check.out`` in the ``analyze/check`` subdirectory, and it will only contains lines in the console output of ``gmx check`` that match "Last frame" and "Coords".  Then, issuing the command:

  .. code-block:: bash

    $ htpolynet analyze -cfg my-analysis.yaml -proj proj-0

  will perform the execution.  Note the ease with which you can reuse this file for all project directories by simply appending project directory names to the ``-proj`` option.  For example, if you have project directories ``proj-0``, ``proj-1``, ``proj-2``, ``proj-3``, and ``proj-4``, you could simply issue the *single* command

  .. code-block:: bash

    $ htpolynet analyze -cfg my-analysis.yaml -proj proj-?

  to run the analysis in each project directory.

2.  ``gmx report-methods``

  ``gmx report-methods`` is a nifty subcommand that generates readable text describing the MD methods used in a given ``tpr`` file.  Suppose you want to use this to remind yourself of the methods used in the densification stage of a particular build, say in project directory ``proj-0``.  You could just issue the ``gmx`` command

  .. code-block:: bash

    $ gmx report-methods -s proj-0/systems/densification/densified-npt.tpr

  Again, very simple.  However, using HTPolyNet, you could make a YAML input file ``another-analysis.yaml``

  .. code-block:: yaml

    - command: report-methods
      subdir: analyze/report-methods
      links:
        - systems/densification/densified-npt.tpr
      options:
        s: densified-npt.tpr
        m: report.tex

  This will generate the tex file ``report.tex`` in the ``analyze/report-methods`` subdirectory of ``proj-0``, if you issue the command

  .. code-block:: bash

    $ htpolynet analyze -cfg another-analysis.yaml -proj proj-0

3. Composites:  You can put as many analysis dictionaries in the single YAML input file.  For example, concatentating the above together would give the single file ``two-analyses.yaml``:

   .. code-block:: yaml

    - command: check
      subdir: analyze/check
      links:
        - systems/densification/densified-npt.trr
      options:
        f: densified-npt.trr
      outfile: check.out
      matchlines: ['Last frame','Coords']
    - command: report-methods
      subdir: analyze/report-methods
      links:
        - systems/densification/densified-npt.tpr
      options:
        s: densified-npt.tpr
        m: report.tex

  And then issuing

  .. code-block:: bash

    $ htpolynet analyze -cfg two-analyses.yaml -proj proj-0

  will run both in series.


These examples should illustrate how powerful ``htpolynet analyze`` is for working with large numbers of parallel systems.

Shortcut analyses
!!!!!!!!!!!!!!!!!

HTPolyNet provides a small number of shortcut analysis dictionaries prepopulated with default values.

1. ``gmx density``

  If an analysis dictionary just contains the lines

  .. code-block:: yaml

    - command: density

  then HTPolyNet assumes the user wants to measure a density profile of all atoms (interative menu option ``0``) along the *z* direction (``-d Z``) in 50 slices (``-sl 50``) by analyzing the trajectory ``equilibrate.trr`` in the subdirectory ``poststim/equilibrate``, generating the free-format output file ``density.xvg`` in the ``analyze/density`` subdirectory.  


2. ``gmx freevolume``

  If an analysis dictionary just contains the lines

  .. code-block:: yaml

    - command: freevolume

  then HTPolyNet assumes the user wants to measure the fractional free volume by analyzing the trajectory ``equilibrate.trr`` and input file ``equilibrate.tpr`` in the subdirectory ``poststim/equilibrate``, generating the free-format output file ``ffv.xvg`` in the ``analyze/freevolume`` subdirectory, and reporting console output matching the following: "Free volume", "Total volume", "Number of molecules", "Density", "Molecular volume Vm assuming homogeneity", "Molecular van der Waals volume assuming homogeneity", and "Fractional free volume".

