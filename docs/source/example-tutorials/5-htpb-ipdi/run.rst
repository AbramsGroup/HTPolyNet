.. _htpb_run:

Running the Build
-----------------

From inside the working directory containing ``5-htpb-ipdi.yaml``:

.. code-block:: console

   $ htpolynet run -diag diagnostics.log 5-htpb-ipdi.yaml &> console.log &

This is by far the longest of the depot examples — plan for **half a
day to overnight** rather than a coffee break.  The stage layout
under ``proj-N/systems/`` is the standard one (``init/``,
``densification/``, ``precure/``, ``iter-K/``, ``capping/``,
``postcure/``, ``final-results/``, plus ``plots/`` and
``profile.json``).

Setup
^^^^^

``htpolynet`` parameterizes the 35 templates discussed in the
:ref:`configuration page <htpb_configuration>`.  Among the
interesting ones:

.. code-block:: text

   INFO> 35 molecules detected in 5-htpb-ipdi.yaml
   INFO>                       explicit: 34
   INFO>     implied by stereochemistry: 1
   INFO>            implied by symmetry: 0
   INFO> OB: generating mol2 from SMILES via RDKit
   INFO> TB: generating mol2 from SMILES via RDKit
   INFO> TBO: generating mol2 from SMILES via RDKit
   INFO> IPD: generating mol2 from SMILES via RDKit
   INFO> AmberTools> generating GAFF parameters from OB.mol2
   ... (4 small constituents parameterized)
   INFO> AmberTools> generating GAFF parameters from A2.mol2
   ... (6 param-stage products parameterized)
   INFO> AmberTools> generating GAFF parameters from A18_I0.mol2
   INFO> AmberTools> generating GAFF parameters from A18_I1.mol2
   ... (16 procession-iteration A18 intermediates parameterized)
   INFO> AmberTools> generating GAFF parameters from DHT.mol2
   INFO> AmberTools> generating GAFF parameters from THT.mol2
   ... (final assembled chains)
   INFO> Generated 35 molecule templates
   INFO> Initial composition is IPD 125, DHT 50, THT 50
   INFO> 100% conversion is 250 bonds

Conformer generation runs next: 6 thermalized conformers of each
``DHT`` and ``THT`` chain at 900 K via short GROMACS NVT MD.  This
takes roughly 5 minutes per conformer (chain × 6 conformers × 2
chains = 12 conformer-generation MD runs), accounting for ~half of
the setup time.

Densification + precure
^^^^^^^^^^^^^^^^^^^^^^^

The 50 densification NPT repeats at 600 K / 10 bar progressively
compact a dilute initial state into a near-melt density of
~0.9-1.0 g/cm³.  Each repeat is 100 ps, and the stage is then extended
further until the density settles -- three extra segments on the run
measured here, ending at 677.2 +/- 0.56 kg/m³.  The full densification
takes ~1h26m of wall clock, which is the single most expensive stage
outside the cure.  Precure adds a 300 ps NPT
preequilibration at 300 K / 1 bar, then an anneal cycle (two
cycles between 300 and 600 K, 200 ps per segment) so the chains can
explore conformational space before cure starts.  Total precure
wall-clock: ~1 hour.

Cure
^^^^

CURE runs until either ``desired_conversion: 0.95`` or
``max_iterations: 150`` is reached.  On the run measured here cure
converges in **9 iterations**.  The per-iteration wall-times are
revealing:

.. list-table::
   :header-rows: 1
   :widths: 15 25 25 35

   * - Iteration
     - Bonds formed
     - Cumulative conversion
     - Wall time
   * - 1
     - 15
     - 0.150
     - 3:26
   * - 2
     - 19
     - 0.340
     - 7:27
   * - 3
     - 10
     - 0.440
     - 6:04
   * - 4
     - 12
     - 0.560
     - 6:52
   * - 5
     - 10
     - 0.660
     - 7:09
   * - 6
     - 10
     - 0.760
     - 12:04
   * - 7
     - 10
     - 0.860
     - 15:23
   * - 8
     - 5
     - 0.910
     - 25:11
   * - 9
     - 4
     - 0.950
     - 23:42

The cure-tail effect is sharp: the first seven iterations together
take ~58 minutes; the last two take another ~49.  By the late iterations
only a handful of hydroxyl / isocyanate pairs are left unbonded,
finding pairs within the bond-search radius requires the
``cure_drag`` step to pull distant atoms together over multiple MD
segments, and each ``cure_drag`` cascade scales with the
inter-atom separation.  ``min_bonds_per_iteration: 10`` is what
keeps the iteration count from blowing up to 50+ at the tail;
raising it further would slightly reduce iteration count but each
iteration would have to drag further-apart atoms together, with
diminishing returns.

Total cure wall-time: ~1h47m.  Capping is trivially fast (0
bonds — all reactive sites that were going to bond did) and runs in
milliseconds.

Postcure
^^^^^^^^

Postcure runs two anneal cycles between 300 K and 600 K (50 ps per
segment) followed by a 200 ps NPT postequilibration at 300 K /
1 bar to let the cured network relax meaningfully before the final
coordinates are written.  Postcure wall-clock: ~7 minutes.

Profile
^^^^^^^

End-of-run stage profile from a representative single-CPU + single-GPU
run:

.. code-block:: text

   Stage                                                   wall      subprocess
   ------------------------------------------------------------------------------
   setup                                                26.94 s          8.05 s
   initialization                                        8.46 s          3.24 s
   densification                                     1h26m07.4s      1h26m03.7s
   precure                                             24m18.1s        24m17.6s
   cure                                               1h47m17.9s            0 ms
     iter-1                                             3m26.4s
     iter-2                                             7m26.9s
     ...
     iter-8                                            25m10.7s
     iter-9                                            23m41.8s
     capping                                              5 ms             0 ms
   postcure                                             6m57.0s         6m56.5s
   final                                                 9.09 s             0 s

Total: **3h45m** on 24 cores.  Of that, gmx-mdrun consumes ~95 % of
the subprocess time; antechamber/parmchk2/tleap account for the rest
of the setup wall.  Wall times scale with core count, so treat these
as a shape rather than a promise: an earlier 16-core run of an
earlier configuration of this example took about 12 hours.

The next page covers the :ref:`results <htpb_results>`.
