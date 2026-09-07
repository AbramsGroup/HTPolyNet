Program Flow
------------

.. mermaid::
   :alt: outer program flow
   :caption: ``htpolynet run`` workflow.

   flowchart TD
       start(["<tt>htpolynet</tt> begins"]):::endpoint
       pre["Make input<br/>molecular structures"]
       setup["Set up reactions and<br/>oligomer templates"]
       topo["Build system topology<br/><i>(e.g., init.top)</i>"]
       coord["Build system coordinates<br/><i>(e.g., init.gro)</i>"]
       md1["<b>MD:</b> Densification and<br/>precure equilibration"]:::md
       cure[["<b>CURE</b>"]]:::cure
       md2["<b>MD:</b> Postcure equilibration"]:::md
       finish(["<tt>htpolynet</tt> ends"]):::endpoint

       start --> pre --> setup --> topo --> coord --> md1 --> cure --> md2 --> finish

       click cure "#cure-section" "Jump to CURE algorithm details"

       classDef endpoint fill:#fff,stroke:#2a7a3f,color:#2a7a3f
       classDef md color:#1f4e9c,stroke:#1f4e9c
       classDef cure fill:#d4edda,stroke:#2a7a3f,color:#1e5128

A basic depiction of the the workflow initiated by ``htpolynet run`` is shown in the figure above.  The first step is generation of the molecular structure data for any monomeric reactants, which ``htpolynet`` does using `RDKit <https://www.rdkit.org/>`_.  Then, based on instructions in the :ref:`configuration file <configuration_files>`, ``htpolynet`` proceeds with setting up all reactions and oligomer templates.  Once these are generated, it then generates the full initial system topology in Gromacs format (that is, it generates a ``top`` file), and an initial set of coordinates (a ``gro`` file).  Since the initial coordinates are built a low density (typically), ``htpolynet`` then performs MD to "densify" the system, followed by any pre-cure equilibration the user would like.  Then the :ref:`CURE algorithm <cure_section>` takes over to generate intermolecular bonds and drive the polymerization.  Once it finishes, ``htpolynet`` conducts any post-cure equilibration the user likes, before saving the final ``top`` and ``gro`` file.  All of this work happens in the project subdirectory of the current directory in which ``run`` is invoked.

This workflow should make clear that the two required tasks of the user are:

1. :ref:`Generating monomer structure files; <molecular_structure_inputs>` and
2. :ref:`Creating an input configuration file. <configuration_files>`

.. _cure_section:

The Connect-Update-Relax-Equilibrate (CURE) algorithm
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. mermaid::
   :alt: block flow diagram of the CURE algorithm
   :caption: Block flow diagram of the CURE algorithm used in ``htpolynet``.

   flowchart TD
       start(["<tt>CURE</tt> begins"]):::beginEnd
       conv{"Conversion<br/>reached?"}
       init["Initialize<br/>search radius"]
       ident["Identify<br/>allowable bonds"]
       found{"Any bonds<br/>found?"}
       inc["Increment<br/>search radius"]
       maxr{"Above max<br/>radius?"}
       longer{"Any longer<br/>than threshold?"}
       drag["Drag"]
       topo1["Update topology"]
       relax1["Relax"]
       equil["Equilibrate"]
       cap["Cap unreacted<br/>groups"]
       topo2["Update topology"]
       relax2["Relax"]
       finish(["<tt>CURE</tt> ends"]):::endEnd

       start --> conv
       conv -- n --> init --> ident --> found
       found -- n --> inc --> maxr
       maxr -- n --> ident
       maxr -- y --> cap
       found -- y --> longer
       longer -- y --> drag --> topo1
       longer -- n --> topo1
       topo1 --> relax1 --> equil --> conv
       conv -- y --> cap
       cap --> topo2 --> relax2 --> finish

       classDef beginEnd fill:#d4edda,stroke:#2a7a3f,color:#1e5128
       classDef endEnd fill:#fadbd8,stroke:#a93226,color:#7b241c

The algorithm used to create new bonds and polymerize a system is called the CURE algorithm, depicted above.  This is just a slightly modified version of a standard search-radius-type algorithm, first used by Li and Strahan to study EPON/DETDA thermosets (:cite:t:`Li2010Crosslinking`).  The CURE algorithm begins by executing a search for new bonds on a frozen system configuration.  Bonds are downselected through a series of filters to arrive at a final set of bonds to form.  If the distance between any pair of "bond-designate" atoms is greater than some threshold (the ``trigger_distance`` parameter in the :ref:`drag subdirective <cure.drag>` of the ``CURE`` directive of a configuration file), a series of MD simulations that slowly bring all to-be-bound atom closer together is performed.  Then the topology is updated, where ``htpolynet`` applies the charges, atom type, and bonded interaction templates from the oligomer template set to each bond.  After the update, a series of relaxation MD simulations bring all bonds to their equilibrium lengths.  Then a short NPT MD simulation equilibrates the overall density before initiating the next CURE iteration.  CURE iterations continue until (a) a desired conversion is reached, or (b) no new allowable bonds are identified.

.. note::

   The drag and relax stages run **unconstrained**, at a 1 fs timestep, because
   the bonds they are manipulating are deliberately far from equilibrium.  Only
   the per-iteration equilibration (and, later, ``postsim``) uses a 2 fs
   timestep with constraints, and as of v2.7.0 it constrains **hydrogen bonds
   only**, with ``lincs_order = 8``.

   Before v2.7.0 it constrained *all* bonds at 2 fs with the GROMACS default
   LINCS accuracy, which is marginal whenever heavy-atom bonds are constrained
   and was fatal for halogenated monomers: a fluorinated bisphenol failed every
   build attempt, dying at the equilibration step several cure iterations in.
   The trap was that the failure surfaces immediately after the relax ladder,
   so the natural diagnosis is that the relax schedule is too coarse -- and it
   is not, since that ladder is unconstrained and cannot be responsible.
   Refining it makes matters worse.  If you see LINCS warnings at
   ``cure_equilibrate``, look at the constraints, not at the ladder.

.. _relax_diagnostics:

What the relax stages report
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As of v2.7.0 the relax ladder reports two diagnostics per CURE iteration.  Both
are pure observation: they change no simulation input and gate nothing.

The per-stage table gains a **Density** column, read from the NPT ``.edr`` each
stage already wrote.  The relax stages are the only above-:math:`T_g`
constant-pressure time in a cure -- roughly 120 ps of it at the defaults,
against a production ladder measured in nanoseconds -- so this is where you can
see whether the box is still densifying when the cure stops.  Drag stages are
excluded, because they run under restraints and their density is not
comparable.

After each ladder, ``htpolynet`` reports **reactive-species mobility**: the rmsd
displacement of atoms that still carry an unused reactive site, and the fraction
of them that moved at least one ``CURE.controls.search_radius`` during the
window.  This is the criterion the fixed-relaxation-window convention was
originally sized against -- that unreacted species diffuse far enough between
reactions to find new partners -- and it has been in use, at 40 ps, since
Varshney and co-workers introduced it in 2008.

That requirement decays over a cure.  Past the gel point the still-reactive
species are bonded into the growing network, and are then topologically
constrained rather than merely slow, so more relaxation time does not restore
their mobility.  A build in which fewer than 25 % of still-reactive atoms cross
a search radius emits a warning, because its later bonds are being chosen from a
nearly frozen neighborhood.  Treat that as a statement about the *validity range
of the protocol*, not as an error: nothing is wrong with the run, but its late
bonds are less well sampled than its early ones.

A run resuming mid-ladder from a checkpoint skips the mobility report, since the
window it could measure is only the tail of the real one.

.. _bondsearch_filters:

Identifying allowable bonds:  Bondsearch filters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A key task of ``htpolynet`` is identifying potential bonds between reactive atoms.  ``htpolynet`` organizes its "bondsearch" in any one CURE iteration along reactions.  Each reaction specified in the :ref:`configuration file <configuration_files>` usually defines one bond by the identies of each reactant and the name of each atom in each reactant.  The bondsearch algorithm begins considering each reaction, and for each one, filtering the set of *possible* bonds using a series of rules, depicted below.

.. figure:: pics/bond_filter.png

   Bondsearch filtering in one CURE iteration.

The outer loop in this figure corresponds to an iteration over reactions.  An iteration consists of downselecting from all atoms to a set of pairs of atoms such that each member corresponds to one of the two atoms referred to in that reaction's bond record.  This selection is also limited to those atoms that still have sacrificial hydrogens.  From this set, any potential bond that would result in a "short-circuit", defined as two monomers that share *more* than one intermonomer bond, is excluded.  Then any potential bond that "pierces" a ring is excluded.  In the case that the polymerization chemistry involves opening of C-C double bonds, any potential bond that results in a "cycle" of C-C bonds of any length is excluded.  This now smaller set of potential bonds are then sorted by length and the topmost potential bonds (shortest) that do not repeat any monomer index are retained.  Then the entire set of potential bonds is considered at once to determine if any cycles of C-C bonds would form, and if so, the longst potential bond in any cycle is disallowed.  Then, for each bond, a random number between 0 and 1 is drawn and compared to its assigned probability; if the random draw is greater than the probability, the bond is disallowed.  Finally, if the current allowable conversion in the iteration is limited by a user directive, the longest bonds beyond this limit are also disallowed.  This final set of bonds is forwarded to the topology update; if this set is of zero length, the topology update immediately hands off to the radius checker.

