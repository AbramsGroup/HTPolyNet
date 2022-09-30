.. _pms_results:

Results
-------

Overall behavior
^^^^^^^^^^^^^^^^

Using ``htpolynet plots`` we can generate a few interesting graphics that help characterize a build.  In this tutorial, we generated a 95%-cure build under ``proj-0``, with diagnostic output in ``diagnostics.log`` and console output to ``console.log``.  

As an exercise, edit ``pMSTY.yaml`` to change the desired cure to 0.50 instead of 0.95.  Then launch a second build:

.. code-block:: bash

    $ htpolynet run -diag diagnostics-low.log pMSTY.yaml &> console-low.log


This will populate the project directory ``proj-1`` (whose name is automatically assigned).  Once it completes, we can generate some plots.

First, we can make plots of the conversion vs. run time and the cure iteration vs. run time:

.. code-block:: console

    $ htpolynet plots -logs diagnostics.log diagnostics-low.log

This generates ``cure-info.png``: 

.. figure:: pics/cure-info.png

    (Left) Conversion vs. wall-clock time; (right) Iteration number vs wall-clock time.

We can see here that the 95\% cure took about 8 and a half minutes of run time (which is not really impressive since this is a **very** small system).  Fully two-thirds of the run time is consumed realizing the final 15\% of the cure.

Second, we can make plots that track the temperature and density throughout the entire build process:

.. code-block:: console

    $ htpolynet plots -proj proj-0 -t p0-traces.png -o p0.csv

This command extracts temperature, density, and potential energy from all Gromacs ``edr`` output files in ``proj-0/`` in the order they were generated, plots them according to a default format, and stores the extracted data in a ``csv`` file.

.. figure:: pics/p0-traces.png 

    (Top) Temperature vs. run time; (middle) Density vs. run time; (bottom) Potential energy vs. run time.  In all panels, vertical lines designate initiations of Gromacs simulations.

From these traces, we can see how little MD time is actually devoted to forming the bonds as compared to relaxing both before and after.  The top plot shows temperature in K vs. time in ps througout the build process.  Vertical lines denote transitions from one step to the next; transitions are very close together in time during the CURE iterations since I'm showing one transition for each drag/relax stage.  The middle plot shows the density trace' note how the density begins at the stipulated low value of 300 kg/m3.  The bottom plot shows the potential energy trace.

In the figure below, we show two renderings of this system.  In each, all bonds between C1 and C2 atoms are shown as grey tubes, and all other bonds are colored by individual unique monomer and made transparent.  On the left is the system just after the precure anneal, where you can see that only **intramolecular** C1 and C2 bonds exist.  On the right is the system after postcure, where you can see chains of -C1-C2- bonds.

.. list-table:: 

    * - .. figure:: pics/hi-pre.png

           Methystyrene liquid before cure.

      - .. figure:: pics/hi.png

           Poly(methyl styrene) after 95% cure.

Details
^^^^^^^

The ``htpolynet run`` invocation in ``run.sh`` runs a high-cure build (95\% conversion) in the ``proj-0`` subdirectory:

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
~~~~~~~~~~~~~~~~~~

.. code-block:: console

    $ cd systems
    $ ls
    capping/         densification/  init/    iter-10/  iter-3/  iter-5/  iter-7/  iter-9/    precure/
    cure_state.yaml  final-results/  iter-1/  iter-2/   iter-4/  iter-6/  iter-8/  postcure/

The ``init/`` directory is where the initial topology and coordinates are generated.  Then in ``densification`` are the files associated with the MD simulations used to densify the initial system.  Next comes the ``precure`` directory, which contains all the results of the precure equilibrations and annealing (if requested).  Next come the iteration directories; here, ten CURE iterations were run.  Then comes the ``capping`` directory where the final topology updates are performed to cap any unreacted monomers (reverting them from their "active" forms to their "proper" forms).  Then comes ``postcure`` equilibration and relaxation.  Finally, in ``final-results`` are the ``top``, ``gro``, and ``grx`` files of the final system; the ``top`` and ``gro`` files can be used right away for Gromacs MD simulations.

``proj-1/plots``
~~~~~~~~~~~~~~~~

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