.. _pms_results:

Results
-------

The first ``htpolynet run`` invocation in ``run.sh`` runs the low-cure build (50\% conversion) in the ``proj-0`` subdirectory.  Let's look at the results in this directory:

.. code-block:: console

    $ cd proj-0
    $ ls
    checkpoint_state.yaml  molecules/  plots/  systems/
    $

The ``yaml`` file is just a checkpoint.  We will consider how to use checkpoints in a separate section.  There are three main subdirectories in all project directories:

* ``molecules/``: This directory has one subdirectory, ``molecules/parameterized`` that contains all files associated with generation and parameterization of all molecules and oligomer templates.  The important files here are those with ``gro``, ``top``, ``itr``, and ``grx`` extensions.

* ``systems/``:  This directory contains all directories where Gromacs runs are conducted.

* ``plots/``: This directory contains some plots generated on the fly.

``proj-0/systems``
^^^^^^^^^^^^^^^^^^

.. code-block:: console

    $ cd systems
    $ ls
    capping/  densification/  final-results/  init/  iter-1/  iter-2/  postcure/  precure/

The ``init/`` directory is where the initial topology and coordinates are generated.  Then in ``densification`` are the files associated with the MD simulations used to densify the initial system.  Next comes the ``precure`` directory, which contains all the results of the precure equilibrations and annealing (if requested).  Next come the iteration directories; here, only two CURE iterations were run (remember the console ouput?) so we see only ``iter-1`` and ``iter-2``.  Then comes the ``capping`` directory where the final topology updates are performed to cap any unreacted monomers (reverting them from their "active" forms to their "proper" forms).  Then comes ``postcure`` equilibration and relaxation.  Finally, in ``final-results`` are the ``top``, ``gro``, and ``grx`` files of the final system.

``proj-0/plots``
^^^^^^^^^^^^^^^^