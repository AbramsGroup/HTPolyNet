.. _configuration_analyze:

Configuration Files for ``htpolynet analyze``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``analyze`` is a subcommand that provides a shortcut interface to ``gmx`` subcommands for performing analyses of MD trajectories generated during ``run`` or ``postsim`` phases.

An ``analyze`` configuration file is structured as a list of dictionaries, each entry of which is a *stage* of the analysis.  Each stage consists of execution of a ``gmx`` analysis and parsing/reporting of its output.  

HTPolyNet recognizes the following analyses:

==============================  ===============  =================
Analysis                        keyword          Description
==============================  ===============  =================
Trajectory info                 ``check``        Runs ``gmx check`` on the indicated trajectory file and reports console output
Local density                   ``density``      Runs ``gmx density`` to measure density profiles
Fractional free volume          ``freevolume``   Runs ``gmx freevolume`` and reports selected console results
==============================  ===============  =================

All analyses 

* Check directives

  A ``check`` directive allows specification of several parameters:

  ================== ===========
  Keyword            Description
  ================== ===========
  ``subdir``         Name of directory in which analysis is conducted relative to the *project* directory
  ``links``          List of files to make symlinks to in the analysis directory.  By default, only ``poststim/equilibrate/equilibrate.trr`` is linked
  ``gromacs``        Gromacs directives; see :ref:`here <gromacs_yaml_directive>`
  ``command``        Gromacs subcommand; by default, this is ``check``, but any other subcommand that needs only 
  ``
  ================== ===========
  

* Density directives
* Freevolume directives

