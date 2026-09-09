.. _ps_run:

Running the Build
-----------------

From inside the working directory containing ``1-polystyrene.yaml``:

.. code-block:: console

   $ htpolynet run -diag diagnostics.log 1-polystyrene.yaml &> console.log &

Tailing ``console.log`` (or ``diagnostics.log``) lets you watch progress
in real time.

Stages
^^^^^^

``htpolynet`` works through the following stages.  The ``setup``,
``initialization``, ``densification``, and ``precure`` stages are the
same as in :ref:`example 0 <liquid_styrene_run>`; the new ones are
``cure``, ``postcure``, and the final-write.

1. **setup.**  ``htpolynet`` reads the SMILES from
   ``constituents.STY.smiles``, generates ``STY.mol2``, parameterizes
   it via AmberTools, and *also* generates and parameterizes all
   chain-expanded dimer/trimer/tetramer templates needed for the cure
   reaction.  All of that is invisible in the directory listing but
   shows up clearly in the profile report at end-of-run.

2. **initialization.**  1000 copies of styrene are placed in a box
   sized for ``initial_density: 300 kg/m³``.

3. **densification.**  Minimization → 10 ps NVT @ 300 K → 200 ps NPT
   @ 300 K, 10 bar, the NPT stage repeating until the density settles
   (one extra segment on a 24-core run; ~4 minutes in total).

4. **precure.**  Pre-equilibration, two 300/600 K annealing cycles,
   post-equilibration.

5. **cure.**  Each iteration writes to ``proj-N/systems/iter-K/``:

   * **bondsearch** identifies eligible C1–C2 pairs within the current
     search radius (grown if needed);
   * **drag** (if configured) pulls candidate pairs together;
   * **update** forms the new bond, deletes sacrificial Hs, and remaps
     the local topology from templates;
   * **relax** + **equilibrate** runs short MD to relieve strain and
     re-equilibrate the box.

   The loop stops when ``desired_conversion`` is reached or
   ``max_iterations`` is exhausted.  After cure, capping converts any
   leftover active STY monomers back to the vinyl form
   (``proj-N/systems/capping/``).

6. **postcure.**  Two more 300/600 K annealing cycles plus a 100 ps NPT.

7. **final.**  ``htpolynet`` writes ``final.gro``, ``final.top``,
   ``final.tpx``, ``final.grx``, and the VMD-friendly
   ``final.viz.psf`` / ``final.viz.tcl`` pair into
   ``proj-N/systems/final-results/``.

Directory layout produced
^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: text

   proj-0/
   ├── plots/
   ├── profile.json           # wall-clock + subprocess timing breakdown
   └── systems/
       ├── init/
       ├── densification/
       ├── precure/
       ├── iter-1/ iter-2/ ...
       ├── capping/
       ├── postcure/
       └── final-results/

Monitoring progress
^^^^^^^^^^^^^^^^^^^

To follow what's happening in real time:

.. code-block:: console

   $ tail -f console.log
   $ # or, just the iteration counters:
   $ grep "^INFO> Iteration" console.log

A clean finish ends with an ``htpolynet runtime ends`` line in the log
and a populated ``final-results/`` directory.  The full run profile
(per-stage wall time plus aggregated subprocess time by tool) is
emitted to the log just before that and dumped to ``profile.json``.
