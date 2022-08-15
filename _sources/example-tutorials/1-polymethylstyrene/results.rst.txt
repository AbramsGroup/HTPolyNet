.. _pms_results:

Results
-------

The first ``htpolynet run`` invocation in ``run.sh`` runs the high-cure build (95\% conversion) in the ``proj-1`` subdirectory.  Let's look at the results in this directory:

.. code-block:: console

    $ cd proj-1
    $ ls
    checkpoint_state.yaml  molecules/  plots/  systems/
    $

The ``yaml`` file is just a checkpoint.  We will consider how to use checkpoints in a separate section.  There are three main subdirectories in all project directories:

* ``molecules/``: This directory has one subdirectory, ``molecules/parameterized`` that contains all files associated with generation and parameterization of all molecules and oligomer templates.  The important files here are those with ``gro``, ``top``, ``itr``, and ``grx`` extensions.

* ``systems/``:  This directory contains all directories where Gromacs runs are conducted.

* ``plots/``: This directory contains some plots generated on the fly.

``proj-1/systems``
^^^^^^^^^^^^^^^^^^

.. code-block:: console

    $ cd systems
    $ ls
    capping/  densification/  final-results/  init/  iter-1/  iter-2/  postcure/  precure/

The ``init/`` directory is where the initial topology and coordinates are generated.  Then in ``densification`` are the files associated with the MD simulations used to densify the initial system.  Next comes the ``precure`` directory, which contains all the results of the precure equilibrations and annealing (if requested).  Next come the iteration directories; here, only two CURE iterations were run (remember the console ouput?) so we see only ``iter-1`` and ``iter-2``.  Then comes the ``capping`` directory where the final topology updates are performed to cap any unreacted monomers (reverting them from their "active" forms to their "proper" forms).  Then comes ``postcure`` equilibration and relaxation.  Finally, in ``final-results`` are the ``top``, ``gro``, and ``grx`` files of the final system.

``proj-0/plots``
^^^^^^^^^^^^^^^^

``HTPolyNet`` generates several plots on the fly during a system build.  

.. code-block:: console

    $ cd ../plots
    $ ls -1
    densification-density.png
    iter-1-cure_equilibrate-density.png
    iter-2-cure_equilibrate-density.png
    iter-3-cure_equilibrate-density.png
    iter-4-cure_equilibrate-density.png
    iter-5-cure_equilibrate-density.png
    iter-6-cure_equilibrate-density.png
    iter-7-cure_equilibrate-density.png
    iter-8-cure_equilibrate-density.png
    iter-9-cure_equilibrate-density.png
    postcure-anneal-T.png
    postcure-postequilibration-density.png
    precure-anneal-T.png
    precure-postequilibration-density.png
    precure-preequilibration-density.png

For example, ``densification-density.png`` indicates that the densification simulation was in fact able to densify the system:

.. figure:: pics/densification-density.png

    ``densification-density.png`` for the high-cure build of polymethylstyrene.

We can check that the annealing cycles were correctly performed from either ``precure-anneal-T.png`` or ``postcure-anneal-T.png``:

.. figure:: pics/postcure-anneal-T.png 

    ``postcure-anneal-T.png``

Finally, we can take a look at the density after the postcure-anneal in ``postcure-postequilibration-density.png``:

.. figure:: pics/postcure-postequilibration-density.png 

    ``postcure-postequilibration-density.png``

Note that the final equilibrated density is about 950 kg/m^3 at 300 K and 1 bar, quite a bit higher than the density of about 800 kb/m^3 liquid styrene at 10 bar and 300 K from the densification simulations.  This result is outside the range expected for `poly(4-methyl styrene) <https://polymerdatabase.com/polymers/poly4-methylstyrene.html>`_ of about 1.01 g/cc, but it's not too suprising given that this is a very small system with a low molecular weight, and it was not very extensively equilibrated. 