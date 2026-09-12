.. _periodic_networks:

Analyzing Trajectories of Periodic Networks
===========================================

A crosslinked network built by ``htpolynet`` is a single covalently connected
molecule whose bonds wrap through the periodic boundaries of the box.  Several
GROMACS tools assume every molecule can be drawn as one contiguous piece, and
for a network that assumption fails.  This page explains the warning you will
see and how to process trajectories for visualization, mean-squared
displacement, and free volume.

Everything below was checked on the DGEBA/PACM system of
:ref:`example 3 <example_tutorials>` with GROMACS 2025.4.

One molecule, on purpose
------------------------

``htpolynet`` writes an entire build as one moleculetype:

.. code-block:: text

   [ molecules ]
   ; Compound #mols
   whole_system 1

A cured network is one molecule, and any unreacted monomers are written into the
same block.  (Topologies from releases before 2.8.0 name it ``None`` instead of
``whole_system``.  The name is cosmetic and nothing depends on it.)

"There were N inconsistent shifts. Check your topology."
--------------------------------------------------------

``gmx trjconv -pbc whole`` builds an unwrapped image of each molecule by walking
its bonds.  For a network that percolates across the periodic boundaries, no
consistent image exists: following a cycle of bonds that crosses the box
returns you to the starting atom shifted by a box vector.  GROMACS reports each
conflict as an inconsistent shift.  ``-pbc mol`` has the same problem, because
it needs whole molecules first.

**The warning is expected and does not indicate a topology error.**  It appears
exactly when the network starts to span the box.  Running
``trjconv -pbc whole`` on the per-iteration equilibration of example 3:

.. list-table::
   :header-rows: 1

   * - Bond conversion
     - Inconsistent shifts
   * - 0.175
     - 0
   * - 0.312
     - 0
   * - 0.510
     - 0
   * - 0.690
     - 2
   * - 0.950
     - 140

The count is zero until the network percolates and then grows with cure.  The
onset between 51% and 69% conversion is consistent with the Flory--Stockmayer
gel point for DGEBA (:math:`f = 2`) with PACM (:math:`f = 4`),
:math:`1/\sqrt{3} \approx 0.58`.

The general rule: **do not try to make a network whole.**  What to do instead
depends on the analysis.

Visualization
-------------

At the end of a build ``htpolynet`` writes ``final.viz.psf`` and
``final.viz.tcl`` alongside the final structure.  The PSF carries the real bond
topology, and the Tcl script hides bonds that cross the periodic boundary, so
the network displays correctly without reconstruction:

.. code-block:: console

   $ vmd final.viz.psf final.gro -e final.viz.tcl

A trajectory can be loaded onto the same PSF.  ``htpolynet make-viz`` regenerates
these files for any ``top``/``gro`` pair; see :doc:`usage`.

Mean-squared displacement
-------------------------

Use ``-pbc nojump``, which makes each atom's trajectory continuous across the
boundaries.  An MSD needs only that; atoms do not have to form whole molecules.
Then **turn off gmx msd's own PBC removal**:

.. code-block:: console

   $ gmx trjconv -pbc nojump -f traj.trr -s topol.tpr -o nojump.xtc
   $ gmx msd -f nojump.xtc -s topol.tpr -normpbc -sel <selection>

.. warning::

   ``gmx msd`` defaults to ``-rmpbc``, which tries to make the network whole in
   every frame.  On a network this does not merely print warnings -- it
   corrupts the result, because atoms get shifted by box vectors between
   frames.  On example 3 the default gave an MSD at 100 ps about 2.4 times
   larger than ``-normpbc`` (0.379 against 0.154 nm²), and the command still
   exits normally.

Two more points.  Per-molecule MSD (``-mol``) is meaningless here, since the
molecule is the whole network; select the atom groups you care about instead.
And avoid ``-nopbc``: in GROMACS 2025.4, ``gmx msd -nopbc`` crashed with a
segmentation fault on this trajectory, while ``-normpbc`` worked.  The packaged
MD parameter files already remove center-of-mass motion
(``comm-mode = Linear``), so no separate correction is needed.

``htpolynet`` does not currently provide an MSD workflow of its own.

Free volume
-----------

``htpolynet analyze`` has a ``freevolume`` stage that runs ``gmx freevolume`` on
the post-build equilibration trajectory; see
:ref:`configuration_analyze`.  It prints the same inconsistent-shift warnings,
but for free volume they turned out to be harmless: on example 3 the fractional
free volume was 0.211 ± 0.007 with default settings, 0.210 ± 0.006 with
``-normpbc``, and 0.211 ± 0.007 on a ``nojump`` trajectory.  Unlike MSD, the
result does not depend on how the network is imaged.

Ignore the per-molecule lines in its output -- ``Number of molecules``,
``Average molar mass`` and ``Molecular volume Vm`` -- which describe the whole
box as a single molecule.

Which trajectory to analyze
---------------------------

Files named like ``4-cure_equilibrate-npt.trr`` are the short equilibrations run
inside every CURE iteration.  The leading ``4`` is the step index, not the
iteration number, and a file with that name exists in each ``iter-N``
directory.  They are cure intermediates, not equilibrated networks.  For
property calculations, use the trajectories from ``htpolynet postsim``, which
runs after the build is finished.
