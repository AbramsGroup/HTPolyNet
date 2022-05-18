The Algorithm
~~~~~~~~~~~~~

HTPolyNet produces Gromacs topology and coordinate files for amorphous, crosslinked polymer systems following the basic steps below.

1. If necessary, GAFF-parameterize any input monomers, and build any intermediate molecules dictated by the polymerization chemistry and GAFF-parameterize those.  Outputs of parameterizations are typically saved in a Library so that repeated parameterizations are not done.
2. Based on a specified composition, generate an initial simulation box, and equilibrate it to a liquid-like density.
3. Perform CURE (Connect-Update-Relax-Equilibrate) iterations to introduce intermonomer bonds according to the input chemistry until a desired conversion is met or an execution threshhold is reached.
4. Finalize by performing any requested post-cure chemistry.
