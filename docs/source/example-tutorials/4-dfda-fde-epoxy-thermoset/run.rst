.. _dfe_run:

Running the Build
-----------------

From inside the working directory containing
``4-dfda-fde-epoxy-thermoset.yaml``:

.. code-block:: console

   $ htpolynet run -diag diagnostics.log 4-dfda-fde-epoxy-thermoset.yaml &> console.log &

Because the chemistry, ``CURE`` controls, ``precure`` cascade, and
``postcure`` cascade are identical to
:ref:`example 3 <pde_run>`, the stage layout and the diagnostic-log
shape are also the same.  Refer back to tutorial 3 for the per-stage
narrative; below are the few items that are specific to this build.

Setup
^^^^^

``setup`` parameterizes the same kind of template set as example 3 —
FDE and DFA primitives, the four primary-to-secondary cure products,
the eight secondary-to-tertiary products, the cap variants, and the
diastereomer pool from FDE's stereocenters.  Furan ring perception adds
a small overhead compared to the all-carbocycle PACM/DGEBA case, but
``antechamber`` and ``tleap`` still complete cleanly thanks to the
atom-mapped SMILES path described in :ref:`monomers <dfe_monomers>`.

Densification, precure, CURE, postcure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These mirror example 3 exactly — same 300 kg/m³ initial density, same
300 ps NPT densification gated on measured density, same 0.5 nm starting search radius, same
0.95 desired conversion (target 380 of 400 bonds), same 0.85 late
threshold, same drag/relax/equilibrate cascades.  Expect a similar
shape: early iterations grab dozens of bonds; the long tail picks up
ones and twos with the search radius growing into the drag regime.

The ``CURE.controls.min_bonds_per_iteration`` knob applies the same
way here as in :ref:`example 3 <pde_run>`.  The default of ``10`` was
chosen using measurements on the DGEBA/PACM system (see the table at
:ref:`pde_run`); the same chemistry runs through this system so the
sweet spot is likely similar.  To explore, drop the knob to ``1`` and
compare ``profile.json`` between runs to see the cost of the
un-batched regime.

Final
^^^^^

After postcure, ``htpolynet`` writes ``final.gro``, ``final.top``,
``final.tpx``, ``final.grx``, the VMD viz pair, and the profile report
into ``proj-N/systems/final-results/``.  The
:ref:`next page <dfe_results>` describes what to look at.
