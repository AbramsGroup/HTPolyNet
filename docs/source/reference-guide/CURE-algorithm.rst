The Connect-Update-Relax-Equilibrate (CURE) Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**This page is under construction.**

Formation of crosslink bonds is a iterative process in HTPolyNet.  We refer to the iteration by the acronym CURE, for "Connect", "Update", "Relax", and "Equilibrate". In detail, the following are the steps of one CURE iteration:

1. **Bond search.**  Bond searching begins with stipulation of a search cutoff radius.  For each reaction, all pairs of reactive atoms that match the reaction pattern are visited to determine if the distance between members of the pair is below cutoff or not.  For each pair below cutoff, a series of secondary filters are applied to determine if the bond is allowable (more on these below).  Every bond that passes the filter is put on a global list.  After all reactions are considered, the resulting list of potential bonds is processed such that only the shortest set of bonds with no repeating residue indices is kept (no residue is allowed to participate in more than one reaction during a single CURE iteration).
2. **Pre-bond directed MD ("dragging").**
3. **Topolgy update.**
4. **Bond relaxation.**
5. **Equilibration.**


Single-pair filters applied in the bond search
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Cutoff**: The atoms in the pair must be closer than the current CURE cutoff distance
2. The bond must not create a second bond from one particular residue to another
3. The bond must not pierce a ring
4. (not yet implemented) The bond must not result in a polyethylene cycle

Filters applied to the set of pairs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. No atom may participate in more than one new bond in any one CURE iteration; if any atom appears more than once in the list of pairs, only the pair with the shortest separation distance is kept.

Dragging
^^^^^^^^

Topology update
^^^^^^^^^^^^^^^

Bond relaxation
^^^^^^^^^^^^^^^

Equilibration
^^^^^^^^^^^^^