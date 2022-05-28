Usage
~~~~~

Installation of the HTPolyNet package gives access to the ``htpolynet`` command.  A typical run can be started like this:

.. code-block:: console

    $ htpolynet run -cfg <my_config.yaml>

where ``<my_config.yaml>`` is replaced by the name of your input configuration file.  A full listing of command-line parameters can be seen using ``htpolynet -h``.

What this does
^^^^^^^^^^^^^^

Invoking ``htpolynet`` always generates a ``proj-n`` subdirectory, where ``n`` is replaced by ``0`` (if the current working directory is empty) or the next integer implied by existing ``proj-0``, ``proj-1``, etc.  This is referred to as a "project" directory.  (If the ``-restart`` flag is appending the the command, it drops into the latest project directory to pick up the most recent checkpoint.)  It then immediately creats ``systems``, ``molecules``, and ``plots`` subdirectories:

* ``molecules/``

  All molecular parameterization files are written to ``molecules/parameterized`` if any are done; otherwise, this directory is populated with any molecule parameterizations (``*.gro``, ``*.top`` / ``*.itp``, and ``*.sea`` files).  (See :ref:`symmetry_equivalence` for a description of ``*.sea`` files.)
* ``systems/``

  * ``init/``: The initial liquid system is built and simulated here.
  * ``iter-n``: The ``n``-th CURE iteration system is simulated here.

* ``plots/``

  A variety of plots are put here.  Currently, only a plot of the density as a function of time from the densification of the initial liquid is put here.  More to come.