.. _pms_run:

Running the Build
=================

Now, in our working directory ``my_pms/2-polymethylstyrene/``, if you haven't already, go ahead and launch via ``./run.sh``.

.. code-block:: console

    $ ./run.sh &

As the ``htpolynet run`` command inside ``run.sh`` indicates, standard output is being redirected from the console to the file ``console.log``.

Console output
^^^^^^^^^^^^^^

In this section, let's go section by section through the console log ``console.log`` to understand what ``HTPolyNet`` is doing.  The first section just shows a banner and indicates where the library is and where the project directory is, then reports on the available Ambertools software, finally ending by stating the name of the configuration file.

.. code-block:: console

    INFO>                                                                    
    INFO>     HTPolyNet 1.0.6                                                
    INFO>     https://abramsgroup.github.io/HTPolyNet/                       
    INFO>                                                                    
    INFO>     Ming Huang                                                     
    INFO>     mh3429@dragons.drexel.edu                                      
    INFO>                                                                    
    INFO>     Cameron F. Abrams                                              
    INFO>     cfa22@drexel.edu                                               
    INFO>                                                                    
    INFO>     Supported in part by Grants W911NF-17-2-0227                   
    INFO>     and W911NF-12-R-0011 from the US Army Research Lab             
    INFO>                                                                    
    INFO> ******************** HTPolyNet runtime begins *********************
    INFO> User library is /home/cfa/htpolynet-tutorials/v1.0.6/2-polymethylstyrene/lib
    INFO> New project in /home/cfa/htpolynet-tutorials/v1.0.6/2-polymethylstyrene/proj-0
    INFO> *************************** Ambertools: ***************************
    INFO> ********************  antechamber (ver. 22.0) *********************
    INFO> ********************        tleap (ver. 22.0) *********************
    INFO> ********************     parmchk2 (ver. 22.0) *********************
    INFO> Configuration: pMSTY.yaml

Next comes the monomer and oligomer template parameterizations.  ``HTPolyNet`` first handles the molecules that are explicitly named in the configuration.  Only ``EMB.mol2`` exists; the other two are generated automatically using their respective ``reaction`` directives.  The three molecules implied by the concept of "chaining" are also generated and parameterized.

.. code-block:: console

    INFO> *********** Templates in proj-0/molecules/parameterized ***********
    INFO> 3 molecules detected in pMSTY.yaml
    INFO>                       explicit: 3    
    INFO>     implied by stereochemistry: 0    
    INFO>            implied by symmetry: 0    
    INFO> AmberTools> generating GAFF parameters from EMB.mol2
    INFO> EMB: 120.19 g/mol
    INFO> AmberTools> generating GAFF parameters from EMB~C1-C2~EMB.mol2
    INFO> EMB~C1-C2~EMB: 238.36 g/mol
    INFO> AmberTools> generating GAFF parameters from EMBCC.mol2
    INFO> EMBCC: 118.17 g/mol
    INFO> 3 molecules implied by chaining
    INFO> AmberTools> generating GAFF parameters from EMB~C1=C2~EMB~C1-C2~EMB.mol2
    INFO> EMB~C1=C2~EMB~C1-C2~EMB: 356.53 g/mol
    INFO> AmberTools> generating GAFF parameters from EMB~C1-C2~EMB~C1=C2~EMB.mol2
    INFO> EMB~C1-C2~EMB~C1=C2~EMB: 356.53 g/mol
    INFO> AmberTools> generating GAFF parameters from EMB~C1-C2~EMB~C1=C2~EMB~C1-C2~EMB.mol2
    INFO> EMB~C1-C2~EMB~C1=C2~EMB~C1-C2~EMB: 474.70 g/mol
    INFO> Generated 6 molecule templates
    INFO> Initial composition is EMB 127
    INFO> 100% conversion is 127 bonds

If we look in ``proj-0/molecules/parameterized`` we'll see the ``gro``, ``itp``, ``top`` and ``grx`` files for each molecule.  (The full parameterizations here were done in the low-cure run since that ``htpolynet`` invocation was first in ``run.sh``.) The first three are Gromacs-specific.  The ``grx`` file contains "extended attributes" of each atom that ``HTPolyNet`` uses internally and are **not** needed for Gromacs.

Next comes initialization of the system topology and coordinates.  Here, using the ``constituents`` directive, ``HTPolyNet`` generates a full system topology and simulation box. The box is filled according to the ``initial_density`` subdirective of the ``densification`` directive.

.. code-block:: console

    INFO> ************** Initialization in proj-0/systems/init **************
    INFO> Topology "init.top" in proj-0/systems/init
    INFO> Initial density: 300.0 kg/m^3
    INFO> Total mass: 2.535e-23 kg
    INFO> Box aspect ratio: 1.0 x 1.0 x 1.0
    INFO> Initial box side lengths: 4.388 nm x 4.388 nm x 4.388 nm
    INFO> Coordinates "init.gro" in proj-0/systems/init
    INFO> Extended attributes "init.grx" in proj-0/systems/init

Next comes a report of the densification of the system.

.. code-block:: console

    INFO> ********** Densification in proj-0/systems/densification **********
    INFO> Running Gromacs: minimization
    INFO> Running Gromacs: nvt ensemble;  10.00 ps,  300.00 K
    INFO> Running Gromacs: npt ensemble; 200.00 ps,  300.00 K,  10.00 bar
    INFO> Current box side lengths: 3.127 nm x 3.127 nm x 3.127 nm
    INFO> Density                      829.27
    INFO> Running-average-Density      726.35
    INFO> Rolling-10-average-Density   828.52
    INFO> Densified coordinates in proj-0/systems/densification/densified-npt.gro

Next comes a report of the precure:

.. code-block:: console

    INFO> **************** Precure in proj-0/systems/precure ****************
    INFO> Running Gromacs: npt ensemble; 200.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 3.133 nm x 3.133 nm x 3.133 nm
    INFO> Density                      824.17
    INFO> Running-average-Density      828.76
    INFO> Rolling-10-average-Density   823.60
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 3.137 nm x 3.137 nm x 3.137 nm
    INFO> Density                      820.71
    INFO> Running-average-Density      830.07
    INFO> Rolling-10-average-Density   829.90

Next we begin the CURE iterations:

.. code-block:: console

    INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> Attempting to form 120 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 1 will generate 39 new bonds
    INFO> Step "cure_relax" initiated on 39 distances (max 0.463 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.405                  0.619
    INFO>          2              0.338                  0.555
    INFO>          3              0.277                  0.504
    INFO>          4              0.230                  0.454
    INFO>          5              0.165                  0.424
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 3.043 nm x 3.043 nm x 3.043 nm
    INFO> Density                      894.71
    INFO> Running-average-Density      843.62
    INFO> Rolling-10-average-Density   892.91
    INFO> Iteration 1 current conversion 0.307 or 39 bonds

This first iteration shows that, with a search radius of 0.5 nm, ``HTPolyNet`` identified 39 allowable bonds.  It then forms them and progresses through the relaxation stages until they are at their correct lengths.  Finally it runs the post-iteration NPT MD equilibration, reporting the resulting box dimensions and density.

Next we proceed through CURE iterations 2 through 13:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 34 new bonds
    INFO> Step "cure_relax" initiated on 34 distances (max 0.493 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.445                  0.660
    INFO>          2              0.386                  0.640
    INFO>          3              0.326                  0.551
    INFO>          4              0.275                  0.509
    INFO>          5              0.214                  0.458
    INFO>          6              0.165                  0.415
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 3.002 nm x 3.002 nm x 3.002 nm
    INFO> Density                      927.98
    INFO> Running-average-Density      889.38
    INFO> Rolling-10-average-Density   927.66
    INFO> Iteration 2 current conversion 0.575 or 73 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 3 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 3 will generate 17 new bonds
    INFO> Step "cure_relax" initiated on 17 distances (max 0.499 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.460                  0.652
    INFO>          2              0.407                  0.631
    INFO>          3              0.329                  0.551
    INFO>          4              0.269                  0.511
    INFO>          5              0.213                  0.463
    INFO>          6              0.167                  0.411
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.975 nm x 2.975 nm x 2.975 nm
    INFO> Density                      950.93
    INFO> Running-average-Density      917.79
    INFO> Rolling-10-average-Density   944.20
    INFO> Iteration 3 current conversion 0.709 or 90 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 4 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 4 will generate 11 new bonds
    INFO> Step "cure_relax" initiated on 11 distances (max 0.465 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.413                  0.622
    INFO>          2              0.349                  0.574
    INFO>          3              0.277                  0.512
    INFO>          4              0.219                  0.461
    INFO>          5              0.156                  0.408
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.995 nm x 2.995 nm x 2.995 nm
    INFO> Density                      930.82
    INFO> Running-average-Density      922.40
    INFO> Rolling-10-average-Density   939.98
    INFO> Iteration 4 current conversion 0.795 or 101 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 5 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 5 will generate 6 new bonds
    INFO> Step "cure_relax" initiated on 6 distances (max 0.495 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.440                  0.660
    INFO>          2              0.379                  0.594
    INFO>          3              0.333                  0.573
    INFO>          4              0.272                  0.503
    INFO>          5              0.212                  0.454
    INFO>          6              0.166                  0.412
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.973 nm x 2.973 nm x 2.973 nm
    INFO> Density                      950.45
    INFO> Running-average-Density      918.67
    INFO> Rolling-10-average-Density   953.45
    INFO> Iteration 5 current conversion 0.843 or 107 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 6 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 6 will generate 3 new bonds
    INFO> Step "cure_relax" initiated on 3 distances (max 0.481 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.417                  0.633
    INFO>          2              0.369                  0.590
    INFO>          3              0.318                  0.540
    INFO>          4              0.271                  0.498
    INFO>          5              0.206                  0.440
    INFO>          6              0.162                  0.399
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.969 nm x 2.969 nm x 2.969 nm
    INFO> Density                      954.32
    INFO> Running-average-Density      930.29
    INFO> Rolling-10-average-Density   951.54
    INFO> Iteration 6 current conversion 0.866 or 110 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 7 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 7 will generate 3 new bonds
    INFO> Step "cure_relax" initiated on 3 distances (max 0.472 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.411                  0.639
    INFO>          2              0.342                  0.588
    INFO>          3              0.282                  0.522
    INFO>          4              0.224                  0.464
    INFO>          5              0.161                  0.414
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.988 nm x 2.988 nm x 2.988 nm
    INFO> Density                      935.63
    INFO> Running-average-Density      892.22
    INFO> Rolling-10-average-Density   941.46
    INFO> Iteration 7 current conversion 0.890 or 113 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 8 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 8 will generate 1 new bond
    INFO> Step "cure_relax" initiated on 1 distance (max 0.349 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.294                  0.526
    INFO>          2              0.251                  0.470
    INFO>          3              0.203                  0.420
    INFO>          4              0.157                  0.400
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.977 nm x 2.977 nm x 2.977 nm
    INFO> Density                      945.79
    INFO> Running-average-Density      922.02
    INFO> Rolling-10-average-Density   942.38
    INFO> Iteration 8 current conversion 0.898 or 114 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 9 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 9 will generate 1 new bond
    INFO> Step "cure_relax" initiated on 1 distance (max 0.459 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.393                  0.627
    INFO>          2              0.342                  0.566
    INFO>          3              0.281                  0.522
    INFO>          4              0.215                  0.443
    INFO>          5              0.160                  0.402
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.990 nm x 2.990 nm x 2.990 nm
    INFO> Density                      933.71
    INFO> Running-average-Density      902.46
    INFO> Rolling-10-average-Density   924.82
    INFO> Iteration 9 current conversion 0.906 or 115 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 10 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Iteration 10 will generate 2 new bonds
    INFO> Step "cure_relax" initiated on 2 distances (max 0.578 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.506                  0.734
    INFO>          2              0.458                  0.690
    INFO>          3              0.401                  0.634
    INFO>          4              0.336                  0.569
    INFO>          5              0.273                  0.503
    INFO>          6              0.216                  0.445
    INFO>          7              0.163                  0.414
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 3.008 nm x 3.008 nm x 3.008 nm
    INFO> Density                      916.97
    INFO> Running-average-Density      897.35
    INFO> Rolling-10-average-Density   928.33
    INFO> Iteration 10 current conversion 0.921 or 117 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 11 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 11 will generate 1 new bond
    INFO> Step "cure_relax" initiated on 1 distance (max 0.423 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.379                  0.616
    INFO>          2              0.314                  0.531
    INFO>          3              0.260                  0.496
    INFO>          4              0.199                  0.437
    INFO>          5              0.153                  0.404
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.990 nm x 2.990 nm x 2.990 nm
    INFO> Density                      933.55
    INFO> Running-average-Density      886.99
    INFO> Rolling-10-average-Density   922.68
    INFO> Iteration 11 current conversion 0.929 or 118 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 12 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 12 will generate 1 new bond
    INFO> Step "cure_relax" initiated on 1 distance (max 0.469 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.408                  0.646
    INFO>          2              0.345                  0.581
    INFO>          3              0.275                  0.503
    INFO>          4              0.218                  0.448
    INFO>          5              0.157                  0.409
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.979 nm x 2.979 nm x 2.979 nm
    INFO> Density                      943.83
    INFO> Running-average-Density      923.67
    INFO> Rolling-10-average-Density   939.60
    INFO> Iteration 12 current conversion 0.937 or 119 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 13 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Iteration 13 will generate 1 new bond
    INFO> Step "cure_relax" initiated on 1 distance (max 0.571 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.505                  0.711
    INFO>          2              0.449                  0.672
    INFO>          3              0.394                  0.622
    INFO>          4              0.328                  0.538
    INFO>          5              0.278                  0.526
    INFO>          6              0.208                  0.451
    INFO>          7              0.155                  0.393
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.976 nm x 2.976 nm x 2.976 nm
    INFO> Density                      946.62
    INFO> Running-average-Density      920.90
    INFO> Rolling-10-average-Density   937.37
    INFO> Iteration 13 current conversion 0.945 or 120 bonds

This meets our desired cure of 95\%, so now ``HTPolyNet`` proceeds to capping, 
and not finding any cappable bonds, proceeds to the postcure:

.. code-block:: console

    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 0 new bonds
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********
    INFO> *************** Postcure in proj-0/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.967 nm x 2.967 nm x 2.967 nm
    INFO> Density                      955.28
    INFO> Running-average-Density      948.42
    INFO> Rolling-10-average-Density   953.16
    INFO> *********** Final data to proj-0/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

This just tells us the final density and where the final results are found. 