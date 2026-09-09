.. _badcy_run:

Running the Build
-----------------

From inside the working directory containing
``6-cyanate-ester.yaml``:

.. code-block:: console

   $ htpolynet run -diag diagnostics.log 6-cyanate-ester.yaml &> console.log &

The stage layout under ``proj-N/systems/`` matches earlier examples
(``init/``, ``densification/``, ``precure/``, ``iter-K/``,
``capping/``, ``postcure/``, ``final-results/``, plus ``plots/`` and
``profile.json`` at the project root) **and adds a new** ``repair/``
**directory** between ``capping/`` and ``postcure/``.  The repair
stage writes its ``repaired.gro``/``repaired.top``/``repaired.tpx``
plus a steepest-descent + short NVT relaxation pair there so the
modified topology has a chance to settle before the postcure MD
ensemble takes over.

Setup
^^^^^

``htpolynet`` parameterizes the 11 templates discussed in the
:ref:`configuration page <badcy_configuration>`:

.. code-block:: text

   INFO> 11 molecules detected in 6-cyanate-ester.yaml
   INFO>                       explicit: 5
   INFO>     implied by stereochemistry: 0
   INFO>            implied by symmetry: 6
   INFO> AmberTools> generating GAFF parameters from BPA.mol2
   INFO> BPA: 228.28 g/mol
   INFO> AmberTools> generating GAFF parameters from TAZ.mol2
   INFO> TAZ: 81.08 g/mol
   INFO> AmberTools> generating GAFF parameters from CYN.mol2
   INFO> CYN: 27.03 g/mol
   INFO> AmberTools> generating GAFF parameters from BPA~O1-C1~TAZ.mol2
   INFO> BPA~O1-C1~TAZ: 307.35 g/mol
   INFO> AmberTools> generating GAFF parameters from BPA~O1-C1~CYN.mol2
   INFO> BPA~O1-C1~CYN: 253.29 g/mol
   ...
   INFO> Generated 11 molecule templates
   INFO> Initial composition is BPA 360, TAZ 240
   INFO> 100% conversion is 720 bonds

The molecular weights are a useful quick sanity check: 228.28 (BPA) +
81.08 (TAZ) – 2 × 1.008 (lost H atoms) = 307.35 (``BPA~O1-C1~TAZ``);
228.28 (BPA) + 27.03 (HCN) – 2 × 1.008 = 253.29
(``BPA~O1-C1~CYN``).

Densification + precure
^^^^^^^^^^^^^^^^^^^^^^^

200 kg/m³ initial density and a 100 ps NPT segment (× 4 repeats) bring
the box to roughly 1.0 g/cm³ before precure.  Precure runs the
preequilibration (200 ps NPT at 300 K, 1 bar) and a short anneal cycle
between 300 and 500 K to relax any high-energy contacts from the
random initial placement.

Cure
^^^^

CURE runs until either ``desired_conversion: 0.90`` or
``max_iterations: 150`` is reached.  On a representative run the cure
converges in nine iterations:

.. list-table::
   :header-rows: 1
   :widths: 20 30 30 20

   * - Iteration
     - Bonds formed
     - Cumulative conversion
     - Wall time
   * - 1
     - 152
     - 0.211
     - 2:18
   * - 2
     - 144
     - 0.411
     - 2:26
   * - 3
     - 133
     - 0.596
     - 2:20
   * - 4
     - 93
     - 0.725
     - 2:10
   * - 5
     - 49
     - 0.793
     - 1:56
   * - 6
     - 31
     - 0.836
     - 1:51
   * - 7
     - 15
     - 0.857
     - 1:44
   * - 8
     - 20
     - 0.885
     - 1:46
   * - 9
     - 12
     - 0.901
     - 2:36

The classic long tail: 80 % cure in 4 iterations, the remaining 10 %
takes another 5.  No capping work because ``etherify`` is the only
cure reaction and the cap stage has nothing to do (no cap directives
in the YAML).

Repair
^^^^^^

After cure converges, the postcure topology-repair stage fires:

.. code-block:: text

   INFO> ************ Postcure repair in proj-0/systems/repair *************
   INFO> triazine_to_cyanate_cap: 57 incomplete TAZ residues identified (171 caps total, 71 free fragments to donate)
   INFO> triazine_to_cyanate_cap: redistributing residual charge -23.0470 across 142 repaired-residue neighbours
   INFO> ******** Postcure repair performed 57 dismantle operations ********
   INFO> Relaxing repaired geometry
   INFO> Running Gromacs: minimization
   INFO> Running Gromacs: nvt ensemble;   5.00 ps,  300.00 K

Decoding the numbers:

* **57 incomplete TAZ residues** out of 240 total — i.e. 183 of the
  240 triazines (~76 %) reached the full 3-bonded state during cure.
  Each incomplete one carries between 0 and 2 bonded BPAs.
* **171 caps total** = 57 × 3.  Each dismantled ring is split into
  three independent -C#N fragments.
* **71 free fragments to donate** = the number of dangling triazine
  C atoms across all incomplete rings.  This is also the number of
  unreacted BPA-OH groups (by atom conservation), so the matching is
  exact and every free fragment finds a home.
* **171 - 71 = 100 in-place caps**: fragments whose ring C atom was
  already bonded to a BPA during cure, so the BPA-O-C bond is
  preserved and only the atom types, bond orders, and angle/dihedral
  parameters update from the templated BPA-O-C#N values.
* **Residual charge ≈ -23 e** distributed across **142 atoms** =
  71 × 2 (the BPA-O atoms newly bonded to free caps, plus the
  CYN-C atoms whose H was deleted).  This is the charge from the
  deleted sacrificial H atoms, redistributed across the heavy-atom
  neighbours so the system stays net-neutral for Ewald.

The repair stage finishes by running a steepest-descent minimization
and a short (5 ps) NVT settle on the modified topology, so any LJ
clashes introduced by physically relocating the free-cap atoms get
relaxed before postcure MD starts.

Postcure
^^^^^^^^

Postcure runs the standard anneal (between 300 K and 500 K, two
cycles) followed by a 100 ps NPT postequilibration at 300 K and 1 bar.
The final density typically lands around 1.1 g/cm³ — a touch lower
than fully cured BADCy (≈ 1.2 g/cm³) because of the residual ``-C#N``
end-groups breaking the network into smaller clusters.

Profile
^^^^^^^

End-of-run stage profile (representative run, 4-core CPU + 1 GPU):

.. code-block:: text

   Stage                                                   wall      subprocess
   ------------------------------------------------------------------------------
   setup                                                 798 ms            0 ms
   initialization                                        6.81 s          2.71 s
   densification                                        5m34.6s         5m34.2s
   precure                                              3m18.7s         3m18.3s
   cure                                                16m51.8s            0 ms
     iter-1                                             2m18s          1m20s
     iter-2                                             2m26s          1m28s
     iter-3                                             2m20s          1m27s
     ...
   repair                                               1m11.1s           11 s
   postcure                                             1m39.3s         1m38.9s
   final                                                 5.07 s            0 s

Total: **28m48s** on 24 cores.  Densification is longer than its
``ps: 200`` suggests because the stage is extended until the density
settles -- two extra segments here, ending at 1113.5 +/- 0.33 kg/m^3.

The cure dominates the run as expected.  The ``repair`` stage's wall
time (~1 min) is split between the surgery itself (~5 s — the
remaining time the minimization and 5-ps NVT settle).  The
``proj-0/profile.json`` file carries the same data in
machine-readable form.

Next is the :ref:`results page <badcy_results>`.
