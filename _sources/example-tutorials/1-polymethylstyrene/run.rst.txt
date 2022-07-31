.. _pms_run:

Running the Build
=================

Now, in our working directory ``my_pms_build``, we are ready to launch HTPolyNet:

.. code-block:: console

    $ cd my_pms_build
    $ ls 
    pMSTY.yaml  lib/
    $ htpolynet run pMSTY.yaml &> info.log &
    [1]
    $

HTPolyNet by default uses the local ``./lib/`` as the molecule library; not including a value for ``-lib`` forces HTPolyNet to use the system library in the ``Library`` subpackage, and if you are making new molecules, they won't be there.  It is also instructive to write console messages to ``info.log``.  Detailed diagnostic messages appear in ``./htpolynet_runtime_diagnostics.log``.  The build can take several minutes, so we are running it in the background.  All the action is happening in ``proj-0`` (and of course being reported on in ``info.log``), so let's look in there.  

Console output
^^^^^^^^^^^^^^

While the build is running, you can monitor its progress using ``tail -f`` on the ``info.log`` file.  A sample of its contents are shown below for an example build, where I show messages related to startup, template parameterization, liquid simulation, the first CURE iteration and the last CURE iteration, and finalization:

.. code-block:: console

    $ more info.log
    INFO> HTPolyNet runtime begins.
    INFO> System library is /home/cfa/Git/HTPolyNet/Library
    INFO> User library is /home/cfa/htpolynet-tests/styrene/my_pms_build/lib
    INFO> Ambertools commands available for HTPolyNet to use:
    INFO>  antechamber (ver.   22.0) at antechamber                                       
    INFO>        tleap (ver.   22.0) at tleap                                             
    INFO>     parmchk2 (ver.   22.0) at parmchk2                                          
    INFO> 
    INFO> Configuration: pMSTY.yaml
    INFO> ********** Generating molecular templates **********
    INFO> ********** 3 molecules explicit in pMSTY.yaml **********
    INFO> EMB, EMB1_1, EMBCC
    INFO> AmberTools is parameterizing EMB.mol2
    INFO> AmberTools is parameterizing EMB1_1.mol2
    INFO> AmberTools is parameterizing EMBCC.mol2
    INFO> ********** 3 molecules implied by chaining **********
    INFO> EMB~C1=C2~EMB1_1, EMB1_1~C1=C2~EMB, EMB1_1~C1=C2~EMB1_1
    INFO> AmberTools is parameterizing EMB~C1=C2~EMB1_1.mol2
    INFO> AmberTools is parameterizing EMB1_1~C1=C2~EMB.mol2
    INFO> AmberTools is parameterizing EMB1_1~C1=C2~EMB1_1.mol2
    INFO> ********** Generated 6 molecule templates **********
    INFO> System initial composition is EMB 100
    INFO> Maximum conversion is 100 bonds.
    INFO> System has 2100 atoms.
    INFO> Initial density: 300.0 kg/m^3
    INFO> Total mass: 1.996e-23 kg
    INFO> Box aspect ratio: 1 x 1 x 1
    INFO> -> Resulting initial box side lengths: 4.052 nm x 4.052 nm x 4.052 nm
    INFO> Generated init.top and init.gro.
    INFO> Conducting initial NPT MD densification simulation of liquid
    INFO> Densified coordinates in npt-1.gro
    INFO> time(ps)                   300.000000
    INFO> density(kg/m^3)            824.783875
    INFO> Running-average-density    772.877910
    INFO> Rolling-average-10         821.042603
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> CURE iteration 1 begins in proj-0/systems/iter-1.
    INFO> Bondsearch using radius 0.5 nm initiated.
    INFO> CURE iteration 1 will generate 23 new bonds.
    INFO> Topology update
    INFO> Relaxation of 23 bonds; max distance 0.476 nm, max 1-4 distance 0.699 nm
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.442                  0.656
    INFO>          2              0.417                  0.632
    INFO>          3              0.376                  0.603
    INFO>          4              0.340                  0.572
    INFO>          5              0.290                  0.514
    INFO>          6              0.259                  0.486
    INFO>          7              0.230                  0.465
    INFO>          8              0.197                  0.447
    INFO>          9              0.167                  0.422
    INFO> Equilibration for 50000 steps at 300 K and 1.0 bar
    INFO>   -> average density 766.016 kg/m^3
    INFO> Current conversion: 0.23 (23/100)
    ... (skipping)
    INFO> CURE iteration 11 begins in proj-0/systems/iter-11.
    INFO> Bondsearch using radius 0.5 nm initiated.
    INFO> Increasing cutoff radius to 0.75 nm
    INFO> CURE iteration 11 will generate 2 new bonds.
    INFO> Prebond dragging initiated on 2 new bonds (max distance 0.617 nm).
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.602
    INFO>          2              0.591
    INFO>          3              0.541
    INFO>          4              0.525
    INFO>          5              0.486
    INFO>          6              0.456
    INFO>          7              0.443
    INFO>          8              0.406
    INFO>          9              0.376
    INFO>         10              0.354
    INFO>         11              0.329
    INFO>         12              0.302
    INFO> Topology update
    INFO> Relaxation of 2 bonds; max distance 0.302 nm, max 1-4 distance 0.562 nm
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.289                  0.547
    INFO>          2              0.261                  0.493
    INFO>          3              0.241                  0.466
    INFO>          4              0.217                  0.459
    INFO>          5              0.184                  0.440
    INFO>          6              0.165                  0.413
    INFO> Equilibration for 50000 steps at 300 K and 1.0 bar
    INFO>   -> average density 812.544 kg/m^3
    INFO> Current conversion: 0.96 (96/100)
    INFO> Current conversion 0.96 exceeds desired conversion 0.95
    INFO> Postcure equilibration.
    INFO> CURE finished.
    INFO> HTPolynet runtime ends.
    
This output should give you some idea of the flow of HTPolyNet.  The first major task is generating molecular templates; in this case, that means the molecules that are named explicitly in the configuration file (either as part of ``initial_composition`` or as a reactant or product of any reaction).  It then builds the simulation box and runs a "densification" MD simulation.  Then, CURE iterations start.  I show here the output from the first iteration and the last iteration only.  Note in the last iteration that the search radius was increased by one increment (0.25 nm), and the new bonds had to be pre-dragged before being formed and relaxed.  The conversion at that point is 96\%, so the algorithm terminates.  Let's look inside the project directory structure to find more results.

Parameterization results
^^^^^^^^^^^^^^^^^^^^^^^^

The first major block of calculations involves parameterizing all required molecular templates.  So after a few minutes, we can inspect the parameterization results while the build continues:

.. code-block:: console

    $ cd proj-0
    $ ls
    molecules/  plots/  systems/
    $ cd molecules/parameterized
    $ ls EMB.*
    EMB.edr  EMB.frcmod  EMB.gro  EMB.grx  EMB.itp  EMB.log  EMB.mol2  EMB.top  EMB.tpr  EMB.trr
    $

What are we seeing here?  These are all the output files for generated by parameterization of the ``EMB`` residue defined in ``lib/molecules/inputs/EMB.mol2``.  six distinct template molecules, along with the EMB monomer.  Each set comprises four files:

1. ``mol2`` -- output of ``antechamber``
2. ``frcmod`` -- output of ``parmchk2``
3. ``gro/top/itp`` -- output of ``parmed`` (``gro`` is also output of Gromacs minimization)
4. ``grx`` -- custom output of HTPolyNet with "extra" atom parameters
5. ``tpr/trr/edr`` -- files associated with a Gromacs-based energy minimization that produces an coordinate-optimized ``gro`` file.

We can see how many molecules have been parameterized (each will have its own ``gro`` file):

.. code-block:: console

    $ ls *.gro 
    'EMB1_1~C1=C2~EMB1_1.gro'  'EMB1_1~C1=C2~EMB.gro'   EMB1_1.gro  'EMB~C1=C2~EMB1_1.gro'   EMBCC.gro   EMB.gro

The filesname in single-quotes are the ones corresponding to the automatically "chain-expanded" reactions. So we have parameterized every reaction product, and this is a sufficient set of templates.

Liquid generation and densification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's back out of the ``molecules/parameterized/`` subdirectory and drop into ``systems/init/``:

.. code-block:: console

    $ cd ../../systems/init
    $ ls
    EMB.gro                 mdout.mdp
    gmx.in                  minimize.mdp
    init.gro                npt-1.cpt
    init.grx                npt-1-density.xvg
    init-minimized.edr      npt-1.edr
    init-minimized.gro      npt-1.gro
    init-minimized.log      npt-1.log
    init-minimized.tpr      npt-1-out.xvg
    init-minimized.trr      npt-1.tpr
    init.top                npt-1.trr
    liquid-densify-npt.mdp
    $

HTPolyNet creates ``init.top`` by merging 100 EMB topologies together and writes it.  ``init.gro`` is created using ``gmx insert-molecules`` using ``EMB.gro`` as input.  ``init.grx`` is created by HTPolyNet to store auxiliary atom attributes.  We can see the results of two ``mdrun`` invocations:

1. ``minimize.mdp`` and ``init-minimized.*`` -- this is just an initial energy minimization of the system; and
2. ``liquid-densify-npt.mdp`` and ``npt-1.*`` -- this is the densification MD simulation that runs until the density is equilibrated at the temperature and pressure specified in the ``mdp`` file.
3. Some ``xvg`` files are created by ``gmx energy``, which is controlled non-interactively using ``gmx.in``.  A plot of density vs. time is generated.

Let's go up out of ``systems/init`` and into ``plots/``:

.. code-block:: console

    $ cd ../../plots
    $ ls
    init-density.png

HTPolyNet automatically generates a plot of density vs. time for the densification simulation:

.. image:: init-density.png

We can see that we've successfully arrived at the liquid-like density of about 810 kg/m\ :sup:`3`. (The actual density is 860, so we are a bit low.)

The first CURE iteration
^^^^^^^^^^^^^^^^^^^^^^^^

The next major part of the build is the first CURE iteration, which is the most expensive because the pair searching is most demanding when the number of reactive atoms is maximal.  So let's back out of ``plots`` and drop into ``systems/iter-1``.  When the iteration completes, there are a _lot_ of files generated.  They can be divided into five "phases" for each CURE iteration:

0. Bond search
1. Dragging
2. Topology update
3. Relaxation
4. Equilibration

Names of files corresponding to phases 0-4 all begin with their respective digits.  Names of files that do not begin with a digit are "auxiliary".  Let's consider the files in the seven states.

Bondsearch files
----------------

Files associated with the bondsearch begin with ``0``:

.. code-block:: console

    $ ls -1 0-*
    0-bondsearch-bonds.csv
    0-bondsearch.gro
    0-bondsearch.grx
    0-bondsearch-input.gro
    0-bondsearch.top

The ``gro`` and ``top`` files are sufficient Gromacs input.  The ``grx`` file contains values of four extra attributes for each atom:

.. code-block:: console

    $ head 0-bondsearch.grx
    globalIdx  z  nreactions reactantName  cycle  cycle_idx  chain  chain_idx
            1  0           0          EMB      0          0     -1         -1
            2  0           0          EMB      0          5     -1         -1
            3  0           0          EMB      0          4     -1         -1
            4  0           0          EMB     -1         -1     -1         -1
            5  0           0          EMB      0          3     -1         -1
            6  0           0          EMB      0          2     -1         -1
            7  0           0          EMB      0          1     -1         -1
            8  1           0          EMB     -1         -1      0          0
            9  1           0          EMB     -1         -1      0          1


``globalIdx`` corresponds to the ``nr`` attribute in the ``[ atoms ]`` directive of a ``top`` file, or the ``atomNum`` attribute of a ``gro`` file; it is just the global atom index.  ``z`` is the current value of the number of available crosslink bonds for that atom.  ``nreactions`` is the number of times the atom has reacted; by default the sum of ``z`` and ``nreactions`` must be a constant.  ``reactantName`` is initialized as the residue name the atom belongs to.  However, as we will see, this attribute is key for communicating which product template maps onto a set of particular residues that react.  ``cycle`` indicates the index of the particular cyclic functional group the atom belongs to, with ``-1`` indicating "none"; here, atoms 1, 2, 3, 5, 6, and 7 all belong to cycle 0.  ``cycle_idx`` is the local index of that atom inside its cycle.  Similarly, ``chain`` indicates the index of the individual chain the atom belongs to; here, atoms 8 and 9 belong to chain 0.  ``chain_idx`` is the local index of that atom inside its chain.  Both ``cycle`` and ``chain`` are maintained as globally unique over the whole system.

The ``csv`` file is a dump of the bonds "DataFrame" showing all the bonds that have been identified and that HTPolyNet intends to implement. If you must, it's good to examine it using ``pandas``:

.. code-block:: console

    $ python
    Python 3.10.5 | packaged by conda-forge | (main, Jun 14 2022, 07:04:59) [GCC 10.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import pandas as pd
    >>> a=pd.read_csv('0-bondsearch-bonds.csv',sep=' ',header=0,index_col=None)
    >>> a.head()
         ai  ri    aj  rj  prob reactantName  order      r       result  allowed  remove-to-uncyclize  lucky  initial_distance
    0  1184  57  1584  76   1.0       EMB1_1      1  0.352  BTRC.passed     True                False   True          0.352211
    1   260  13   933  45   1.0       EMB1_1      1  0.354  BTRC.passed     True                False   True          0.353531
    2   491  24  1290  62   1.0       EMB1_1      1  0.363  BTRC.passed     True                False   True          0.363341
    3  1499  72   996  48   1.0       EMB1_1      1  0.369  BTRC.passed     True                False   True          0.368711
    4  2003  96   681  33   1.0       EMB1_1      1  0.371  BTRC.passed     True                False   True          0.371210
    >>>

``ai`` and ``aj`` are the global atom indices for each bond-designate; ``ri`` and ``rj`` are their respective residue indices.  For example, the first bond will join residue 57 to 76.  ``prob`` is the a priori probability that the bond was permitted, and it is just whatever is in the ``probability`` field of the reaction that defines this bond in the input file.  ``reactantName`` indicates the the product template of the bond.  ``order`` is the bond order, again just from the reaction input.  ``r`` is the instantaneous interatomic distance in nm.  ``result`` is a code indicating whether or not it passed the bond filters (all bonds that are saved to this file survived, so they will all have ``BTRC.passed`` as a ``result``.)  ``allowed``, ``remove-to-uncyclize``, and ``lucky`` are Booleans used in the bond filtering process, and their particular meanings are not crucial to understand here.  ``initial_distance`` is the same as ``r`` but it is the result of computing distances by a different module and is retained for consistency-checking.

Dragging files
--------------

Files associated with prebond dragging begin with ``1``.  However, because no bond-designate length exceeded 0.5 nm, no dragging is triggered.  So the build proceeds to the topology update.

Topology update files
---------------------

Files associated with the topology update process begin with a ``2``:

.. code-block:: console

    $ ls -1 2*
    2-update-complete-bonds.csv
    2-update-complete.gro
    2-update-complete.grx
    2-update-complete.top
    2-update-idx-mapper.dat
    2-update-resid-graph.json

All files here represent **outputs** of the topology update.  Let's look at the ``2-update-idx-mapper.dat``:

.. code-block:: console

    $ tail 2-update-idx-mapper.dat 
    2089 1977
    2090 1978
    2091 1979
    2092 1980
    2093 1981
    2094 1982
    2095 1983
    2096 1984
    2098 1985
    2099 1986

The purpose of this file is very simple:  The first column are atom indices **before** topology update, and the second column are indices **after** topology update.  Remember that topology updating deletes sacrificial hydrogens, which means atoms are reindexed (since Gromacs requires sequential atom indexes).  This file allows us to match any atoms in pre-update ``gro`` and ``top`` files to those that exist downstream of a topology update.  Note that I've chosen to show a ``tail`` of this file to highlight the largest index differences.  The post-update indexes also appear in the ``csv`` file showing all bonds.

The file ``2-update-complete-bonds.csv`` is just the initial ``0-bondsearch-bonds.csv``, except all the atom indexes have been updated according to the index mapper described above.  The last line in this file reports the new bond with the longest initial length:

.. code-block:: console

    $ python
    Python 3.10.5 | packaged by conda-forge | (main, Jun 14 2022, 07:04:59) [GCC 10.3.0] on linux
    Type "help", "copyright", "credits" or "license" for more information.
    >>> import pandas as pd
    >>> a=pd.read_csv('2-update-complete-bonds.csv',sep=' ',header=0,index_col=None)
    >>> a.head()
         ai  ri    aj  rj  prob reactantName  order      r       result  allowed  remove-to-uncyclize  lucky  initial_distance
    0  1157  57  1547  76   1.0       EMB1_1      1  0.352  BTRC.passed     True                False   True          0.351936
    1   254  13   911  45   1.0       EMB1_1      1  0.354  BTRC.passed     True                False   True          0.353531
    2   478  24  1258  62   1.0       EMB1_1      1  0.363  BTRC.passed     True                False   True          0.363341
    3  1464  72   972  48   1.0       EMB1_1      1  0.369  BTRC.passed     True                False   True          0.368711
    4  1959  96   665  33   1.0       EMB1_1      1  0.371  BTRC.passed     True                False   True          0.371199
    >>> a.tail()
          ai  ri    aj  rj  prob reactantName  order      r       result  allowed  remove-to-uncyclize  lucky  initial_distance
    18  1298  64  1238  61   1.0       EMB1_1      1  0.429  BTRC.passed     True                False   True          0.428661
    19   334  17   645  32   1.0       EMB1_1      1  0.439  BTRC.passed     True                False   True          0.438782
    20   991  49   768  38   1.0       EMB1_1      1  0.440  BTRC.passed     True                False   True          0.439869
    21   294  15  1117  55   1.0       EMB1_1      1  0.470  BTRC.passed     True                False   True          0.470325
    22   747  37   275  14   1.0       EMB1_1      1  0.476  BTRC.passed     True                False   True          0.475644
    >>>

The bonds in ``2-update-complete-bonds.csv`` are the same as those in ``0-bondsearch-bonds.csv`` except with updated atom indices.  Note for instance that the first bond still indicates a linkage between residues 57 and 76.

Again, the ``gro`` and ``top`` are proper Gromacs inputs, and the ``grx`` file tabulates all the "extended" attributes first presented when describing ``0-bondsearch.grx``.  The ``json`` file represents the graph structure of the network on a resid basis in JSON format.

Relaxation files
----------------

Files that begin with a ``3`` correspond to bond relaxation stages.  In this example, six stages are run by virture of the bond-designate with the longest bond length (0.499 nm) and the ``relax_increment`` of 0.075 nm.  Each stage produces 22 output files: the bonds ``csv``, the ``gro`` / ``grx`` / ``top`` that initializes the first stage, and then the 17 **outputs** from the minimization (5), nvt (6), and npt (6) sub-stages:

.. code-block:: console

    $ ls 3-*
    3-relax-stage-1-bonds.csv
    3-relax-stage-1.gro
    3-relax-stage-1.grx
    3-relax-stage-1-min.edr
    3-relax-stage-1-min.gro
    3-relax-stage-1-min.log
    3-relax-stage-1-min.tpr
    3-relax-stage-1-min.trr
    3-relax-stage-1-npt.cpt
    3-relax-stage-1-npt.edr
    3-relax-stage-1-npt.gro
    3-relax-stage-1-npt.log
    3-relax-stage-1-npt.tpr
    3-relax-stage-1-npt.trr
    3-relax-stage-1-nvt.cpt
    3-relax-stage-1-nvt.edr
    3-relax-stage-1-nvt.gro
    3-relax-stage-1-nvt.log
    3-relax-stage-1-nvt.tpr
    3-relax-stage-1-nvt.trr
    3-relax-stage-1.top
    ...
    3-relax-stage-6-bonds.csv
    3-relax-stage-6.gro
    3-relax-stage-6.grx
    3-relax-stage-6-min.edr
    3-relax-stage-6-min.gro
    3-relax-stage-6-min.log
    3-relax-stage-6-min.tpr
    3-relax-stage-6-min.trr
    3-relax-stage-6-npt.cpt
    3-relax-stage-6-npt.edr
    3-relax-stage-6-npt.gro
    3-relax-stage-6-npt.log
    3-relax-stage-6-npt.tpr
    3-relax-stage-6-npt.trr
    3-relax-stage-6-nvt.cpt
    3-relax-stage-6-nvt.edr
    3-relax-stage-6-nvt.gro
    3-relax-stage-6-nvt.log
    3-relax-stage-6-nvt.tpr
    3-relax-stage-6-nvt.trr
    3-relax-stage-6.top

The attenuation is managed by the sequential ``top`` files.  Let's look at the entry for a particular bond (between atoms 581 and 1033) in each stage's ``top`` file's ``[ bonds ]`` directive:

.. code-block:: console

    $ grep "^275 747" 3-relax-stage-?.top|awk '{if ($3==1) print $0}'
    3-relax-stage-1.top:275 747 1 0.4398777993494193 27977.013333333332
    3-relax-stage-2.top:275 747 1 0.4041118244307419 55954.026666666665
    3-relax-stage-3.top:275 747 1 0.36834584951206445 83931.04
    3-relax-stage-4.top:275 747 1 0.33257987459338706 111908.05333333333
    3-relax-stage-5.top:275 747 1 0.2968138996747096 139885.06666666668
    3-relax-stage-6.top:275 747 1 0.2610479247560322 167862.08
    3-relax-stage-7.top:275 747 1 0.22528194983735483 195839.09333333332
    3-relax-stage-8.top:275 747 1 0.18951597491867744 223816.10666666666
    3-relax-stage-9.top:275 747 1 0.15375 251793.12
    $

In a ``[ bonds ]`` topology directive, the 4th and 5th columns are THE ``b0`` and ``kt`` harmonic bond parameters.  In the stage-9 ``top``, we see these parameters at their proper force-field values for a C-C single bond.  Notice how the value of the distance parameter ``b0`` begins at a large initial value and linearly decreases toward the target (but never by *more* than an increment of 0.075 nm), while the spring constant ``kt`` starts low and increases linearly toward its target.  

Equilibration files
-------------------

Files associated with final equilibration of the bonded system at the end of one CURE iteration begin with a ``4``:

.. code-block:: console

    $ ls 4-*
    4-equilibrate-bonds.csv
    4-equilibrate-complete-bonds.csv
    4-equilibrate-complete.gro
    4-equilibrate-complete.grx
    4-equilibrate-complete.top
    4-equilibrate.gro
    4-equilibrate.grx
    4-equilibrate.mdp
    4-equilibrate-post.cpt
    4-equilibrate-post.edr
    4-equilibrate-post.gro
    4-equilibrate-post.log
    4-equilibrate-post.tpr
    4-equilibrate-post.trr
    4-equilibrate.top

Files with the simple prefix ``4-equilibrate`` represent inputs to the Gromacs run.  Files with the prefixs ``4-equilibrate-post`` are the raw Gromacs mdrun outputs, and the files with the prefix ``4-equilibrate-complete`` represent the Gromacs outputs read back in to HTPolyNet and processed.  This set of ``complete`` files are copied to the next CURE iteration directory as the set of ``0-connect`` files.

Subsequent CURE iterations
^^^^^^^^^^^^^^^^^^^^^^^^^^

The primary result of a CURE iteration is the calculated conversion, or the fraction of the maximum number of crosslink bonds possible, based on the initial composition and reaction stoichiometries, that have formed up to that point. If this fraction is below the value associated with the ``CURE_desired_conversion`` option, then a new iteration is begun.  This involves creating the next ``iter-n/`` directory under ``systems/``, and copying over the prior iteration's ``4-equilibrate-complete.top/gro/grx`` files onto the new ``0-bondsearch.top/gro/grx`` files.  At the beginning of any CURE iterations, the maximum number of new bonds required to reach the desired conversion is calculated and used as a limit in creating new bonds, so that the desired conversion is hit exactly.  

The number of CURE iterations needed to reach a specified conversion is never deterministic because of the randomness inherent in the inter-stage and post-bonding MD simulations.  In this particular instance, a total of 9 CURE iterations were required to reach 0.95 conversion.  Files for each iteration's directory follow the same naming convention explained for the first iteration.

Post-cure reactions, equilibration, and finalization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After iteration 9, when the conversion specification is satisfied, HTPolyNet progresses to the post-cure stage.  The directory ``systems/postcure`` is created and the final outputs from the last CURE iterations are copied here.  If there were any monomers that had not yet reacted (here there are not), then the EMBCC reaction would be used to revert them back to double bonds, followed by an equilibration. After the equilibration, HTPolyNet generates the final files ``7-final.top/gro/grx``. 
    
Overall behavior
^^^^^^^^^^^^^^^^

If the build is run with ``--loglevel debug`` indicated on the command-line, the log file will contain a lot of information that can be used to characterize the efficiency of the build process.  The ``HTPolyNet.plot`` module has a method ``cure_graph`` that can be used to generate plots showing the conversion vs. run time in hours, and the iteration number vs. run time in hours.  Generating this plot from the directory the log file is in can be done using an interactive python session:

.. code-block:: python

    >>> from HTPolyNet.plot import cure_graph
    >>> cure_graph(['my_build.log'],xmax=20.)

We ran 10 independent system builds of 100 monomers each using the provided ``mol2`` and ``yaml`` input files; they generated the logs ``0.log``, ``1.log``, ..., ``9.log``.  The plot below was made using:

.. code-block:: python

    >>> from glob import glob
    >>> from HTPolyNet.plot import cure_graph
    >>> cure_graph([glob('[0-9].log')],xmax=0.3)

.. image:: iter-graph.png

In this case, on a moderately slow workstation, these builds took 10-15 minutes to reach 0.95 conversion, usually in 9 iterations.

Below is a trace of the density vs time as a concatenation of the sequence of all NPT MD simulations, beginning with the initial densification, passing through all drag/relaxationg/equilibrations in each iteration, and concluding with the final equilibration:

.. image:: all-density.png

It is clear that during the post-bond relaxations, density drops to 700 kg/m3, but this is because the post-bond relaxations are all run at 600 K.  The equilibrations at 300 K all bring the system back to approx. 900 kg/m3.

