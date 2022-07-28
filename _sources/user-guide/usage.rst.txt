Usage
~~~~~

Installation of the HTPolyNet package gives access to the ``htpolynet`` command.  A typical run can be started like this:

.. code-block:: console

    $ htpolynet run <my_config.yaml>

where ``<my_config.yaml>`` is replaced by the name of your input configuration file.  A full listing of command-line parameters can be seen using ``htpolynet -h``::

  usage: htpolynet [-h] [-lib LIB] [-log LOG] [-restart] [--force-parameterization] [--force-checkin]
                 [--loglevel LOGLEVEL]
                 command config

positional arguments:
  command               command (info, parameterize, run)
  config                input configuration file in YAML format

options:
  -h, --help            show this help message and exit
  -lib LIB              local user library of molecular structures and parameterizations
  -log LOG              diagnostic log file
  -restart              restart in latest proj dir
  --force-parameterization
                        force GAFF parameterization of any input mol2 structures
  --force-checkin       force check-in of any generated parameter files to the system library
  --loglevel LOGLEVEL   Log level for messages written to diagnostic log (debug|info)

What this does
^^^^^^^^^^^^^^

Invoking ``htpolynet`` always generates a ``proj-n`` subdirectory, where ``n`` is replaced by ``0`` (if the current working directory is empty) or the next integer implied by existing ``proj-0``, ``proj-1``, etc.  This is referred to as a "project" directory.  (If the ``-restart`` flag is appending the the command, it drops into the latest project directory to pick up the most recent checkpoint.)  Under the project directory, it then immediately creates ``systems``, ``molecules``, and ``plots`` subdirectories:

* ``molecules/``

  All molecular parameterization files are written to ``molecules/parameterized`` if any are done; otherwise, this directory is populated with any molecule parameterizations (``*.gro``, ``*.top`` / ``*.itp``, and ``*.sea`` files).  (See :ref:`symmetry_equivalence` for a description of ``*.sea`` files.)
* ``systems/``

  * ``init/``: The initial liquid system is built and simulated here.
  * ``iter-n``: The ``n``-th CURE iteration system is simulated here.
  * ``final``: The final output coordinate and topology files are here.

* ``plots/``

  A variety of plots are put here.  Currently, only a plot of the density as a function of time from the densification of the initial liquid is put here.

HTPolyNet reports on the progress of execution by messages to the terminal console.  It also generates a lot of diagnostic output, which by default is saved in the file ``htpolynet_runtime_diagnostics.log``.  The diagnostic output is crucial for us to understand why an execution fails, so it is a good idea to keep it.  Examples of the console output and diagnostic logs are discussed in the :ref:`example_tutorials`.