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
    INFO>     HTPolyNet 1.0.5                                                
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
    INFO> User library is /home/cfa/htpolynet-tests/v1.0.5/2-polymethylstyrene/lib
    INFO> New project in /home/cfa/htpolynet-tests/v1.0.5/2-polymethylstyrene/proj-0
    INFO> *************************** Ambertools: ***************************
    INFO>   antechamber (ver. 22.0) at antechamber                                        
    INFO>         tleap (ver. 22.0) at tleap                                              
    INFO>      parmchk2 (ver. 22.0) at parmchk2                                           
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
    INFO> Initial composition is EMB 100
    INFO> 100% conversion is 100 bonds

If we look in ``proj-0/molecules/parameterized`` we'll see the ``gro``, ``itp``, ``top`` and ``grx`` files for each molecule.  (The full parameterizations here were done in the low-cure run since that ``htpolynet`` invocation was first in ``run.sh``.) The first three are Gromacs-specific.  The ``grx`` file contains "extended attributes" of each atom that ``HTPolyNet`` uses internally and are **not** needed for Gromacs.

Next comes initialization of the system topology and coordinates.  Here, using the ``constituents`` directive, ``HTPolyNet`` generates a full system topology and simulation box. The box is filled according to the ``initial_density`` subdirective of the ``densification`` directive.

.. code-block:: console

    INFO> ************** Initialization in proj-0/systems/init **************
    INFO> Topology "init.top" in proj-0/systems/init
    INFO> Initial density: 300.0 kg/m^3
    INFO> Total mass: 1.996e-23 kg
    INFO> Box aspect ratio: 1.0 x 1.0 x 1.0
    INFO> Initial box side lengths: 4.052 nm x 4.052 nm x 4.052 nm
    INFO> Coordinates "init.gro" in proj-0/systems/init
    INFO> Extended attributes "init.grx" in proj-0/systems/init

Next comes a report of the densification of the system.

.. code-block:: console

    INFO> ********** Densification in proj-0/systems/densification **********
    INFO> Running Gromacs: minimization
    INFO> Running Gromacs: nvt ensemble;  10.00 ps,  300.00 K
    INFO> Running Gromacs: npt ensemble; 200.00 ps,  300.00 K,  10.00 bar
    INFO> Current box side lengths: 2.885 nm x 2.885 nm x 2.885 nm
    INFO> Density                      830.83
    INFO> Running-average-Density      739.88
    INFO> Rolling-10-average-Density   833.65
    INFO> Densified coordinates in proj-0/systems/densification/densified-npt.gro

Next comes a report of the precure:

.. code-block:: console

    INFO> **************** Precure in proj-0/systems/precure ****************
    INFO> Running Gromacs: npt ensemble; 200.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.896 nm x 2.896 nm x 2.896 nm
    INFO> Density                      821.96
    INFO> Running-average-Density      827.60
    INFO> Rolling-10-average-Density   823.81
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.880 nm x 2.880 nm x 2.880 nm
    INFO> Density                      835.61
    INFO> Running-average-Density      830.84
    INFO> Rolling-10-average-Density   829.83

Next we begin the CURE iterations:

.. code-block:: console

    INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> Attempting to form 95 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 1 will generate 28 new bonds
    INFO> Step "cure_relax" initiated on 28 distances (max 0.452 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.405                  0.631
    INFO>          2              0.341                  0.566
    INFO>          3              0.277                  0.500
    INFO>          4              0.214                  0.454
    INFO>          5              0.164                  0.418
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.819 nm x 2.819 nm x 2.819 nm
    INFO> Density                      886.28
    INFO> Running-average-Density      839.07
    INFO> Rolling-10-average-Density   881.52
    INFO> Iteration 1 current conversion 0.280 or 28 bonds

This first iteration shows that, with a search radius of 0.5 nm, ``HTPolyNet`` identified 28 allowable bonds.  It then forms them and progresses through the relaxation stages until they are at their correct lengths.  Finally it runs the post-iteration NPT MD equilibration, reporting the resulting box dimensions and density.

Next we proceed through CURE iterations 2 through 9:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 24 new bonds
    INFO> Step "cure_relax" initiated on 24 distances (max 0.495 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.444                  0.654
    INFO>          2              0.401                  0.618
    INFO>          3              0.319                  0.558
    INFO>          4              0.274                  0.516
    INFO>          5              0.217                  0.456
    INFO>          6              0.162                  0.414
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.766 nm x 2.766 nm x 2.766 nm
    INFO> Density                      934.62
    INFO> Running-average-Density      877.79
    INFO> Rolling-10-average-Density   930.05
    INFO> Iteration 2 current conversion 0.520 or 52 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 3 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 3 will generate 17 new bonds
    INFO> Step "cure_relax" initiated on 17 distances (max 0.476 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.411                  0.638
    INFO>          2              0.362                  0.592
    INFO>          3              0.286                  0.508
    INFO>          4              0.225                  0.461
    INFO>          5              0.170                  0.414
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.754 nm x 2.754 nm x 2.754 nm
    INFO> Density                      944.87
    INFO> Running-average-Density      912.50
    INFO> Rolling-10-average-Density   941.96
    INFO> Iteration 3 current conversion 0.690 or 69 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 4 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 4 will generate 13 new bonds
    INFO> Step "cure_relax" initiated on 13 distances (max 0.497 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.455                  0.679
    INFO>          2              0.382                  0.610
    INFO>          3              0.329                  0.555
    INFO>          4              0.268                  0.512
    INFO>          5              0.217                  0.455
    INFO>          6              0.163                  0.414
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.756 nm x 2.756 nm x 2.756 nm
    INFO> Density                      940.21
    INFO> Running-average-Density      916.85
    INFO> Rolling-10-average-Density   944.55
    INFO> Iteration 4 current conversion 0.820 or 82 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 5 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 5 will generate 4 new bonds
    INFO> Step "cure_relax" initiated on 4 distances (max 0.499 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.436                  0.648
    INFO>          2              0.393                  0.593
    INFO>          3              0.332                  0.560
    INFO>          4              0.267                  0.510
    INFO>          5              0.208                  0.447
    INFO>          6              0.163                  0.414
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.770 nm x 2.770 nm x 2.770 nm
    INFO> Density                      925.04
    INFO> Running-average-Density      907.14
    INFO> Rolling-10-average-Density   926.97
    INFO> Iteration 5 current conversion 0.860 or 86 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 6 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 6 will generate 3 new bonds
    INFO> Step "cure_relax" initiated on 3 distances (max 0.481 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.430                  0.660
    INFO>          2              0.375                  0.613
    INFO>          3              0.323                  0.547
    INFO>          4              0.263                  0.477
    INFO>          5              0.209                  0.460
    INFO>          6              0.156                  0.402
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.771 nm x 2.771 nm x 2.771 nm
    INFO> Density                      923.83
    INFO> Running-average-Density      900.40
    INFO> Rolling-10-average-Density   920.28
    INFO> Iteration 6 current conversion 0.890 or 89 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 7 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 7 will generate 1 new bond
    INFO> Step "cure_relax" initiated on 1 distance (max 0.406 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.337                  0.579
    INFO>          2              0.301                  0.528
    INFO>          3              0.259                  0.493
    INFO>          4              0.207                  0.445
    INFO>          5              0.160                  0.411
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.767 nm x 2.767 nm x 2.767 nm
    INFO> Density                      927.42
    INFO> Running-average-Density      889.26
    INFO> Rolling-10-average-Density   920.36
    INFO> Iteration 7 current conversion 0.900 or 90 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 8 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Iteration 8 will generate 3 new bonds
    INFO> Step "cure_drag" initiated on 3 distances (max 0.688 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.651
    INFO>          2              0.587
    INFO>          3              0.541
    INFO>          4              0.504
    INFO>          5              0.444
    INFO>          6              0.400
    INFO>          7              0.349
    INFO>          8              0.305
    INFO> Step "cure_relax" initiated on 3 distances (max 0.305 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.269                  0.500
    INFO>          2              0.217                  0.448
    INFO>          3              0.166                  0.413
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.759 nm x 2.759 nm x 2.759 nm
    INFO> Density                      935.93
    INFO> Running-average-Density      877.84
    INFO> Rolling-10-average-Density   936.27
    INFO> Iteration 8 current conversion 0.930 or 93 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 9 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Iteration 9 will generate 1 new bond
    INFO> Step "cure_drag" initiated on 1 distance (max 0.712 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.657
    INFO>          2              0.608
    INFO>          3              0.562
    INFO>          4              0.503
    INFO>          5              0.454
    INFO>          6              0.415
    INFO>          7              0.351
    INFO>          8              0.301
    INFO> Step "cure_relax" initiated on 1 distance (max 0.301 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.263                  0.493
    INFO>          2              0.215                  0.451
    INFO>          3              0.153                  0.406
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.753 nm x 2.753 nm x 2.753 nm
    INFO> Density                      941.30
    INFO> Running-average-Density      904.43
    INFO> Rolling-10-average-Density   940.74
    INFO> Iteration 9 current conversion 0.940 or 94 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 10 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Radius increased to 1.0 nm
    INFO> Iteration 10 will generate 2 new bonds
    INFO> Step "cure_drag" initiated on 2 distances (max 0.938 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.902
    INFO>          2              0.818
    INFO>          3              0.768
    INFO>          4              0.701
    INFO>          5              0.644
    INFO>          6              0.589
    INFO>          7              0.540
    INFO>          8              0.472
    INFO>          9              0.412
    INFO>         10              0.360
    INFO>         11              0.306
    INFO> Step "cure_relax" initiated on 2 distances (max 0.306 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.253                  0.478
    INFO>          2              0.207                  0.428
    INFO>          3              0.158                  0.397
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 2.802 nm x 2.802 nm x 2.802 nm
    INFO> Density                      892.66
    INFO> Running-average-Density      865.98
    INFO> Rolling-10-average-Density   893.69
    INFO> Iteration 10 current conversion 0.960 or 96 bonds

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
    INFO> Current box side lengths: 2.731 nm x 2.731 nm x 2.731 nm
    INFO> Density                      964.01
    INFO> Running-average-Density      949.80
    INFO> Rolling-10-average-Density   952.64
    INFO> *********** Final data to proj-0/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

This just tells us the final density and where the final results are found.  If we look there, we see:

.. code-block:: console

    $ ls -l proj-0/systems/final-results
    final.gro  final.grx  final.top

Now, with the ``gro`` and ``top`` file, you can run whatever Gromacs simulation you like with this system.

