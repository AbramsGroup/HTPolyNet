.. _ve_run:

Running the Build
=================

Now, in our working directory ``my_vinyl_ester/3-bisgma-styrene-thermoset/``, if you haven't already, go ahead and launch via ``./run.sh``.

.. code-block:: console

    $ ./run.sh

As the two ``htpolynet run`` commands indicate, standard output is being redirected from the console to the file ``lo.log`` for the low-cure run and ``hi.log`` for the high-cure run.  Let's consider the contents of ``lo.log``.

First the banner and the runtime initialization messages:

.. code-block:: console

    INFO>                                                                    
    INFO>     HTPolyNet 0.0.1                                                
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
    INFO> User library is /home/cfa/htpolynet-tests/testset/3-bisgma-styrene-thermoset/lib
    INFO> New project in /home/cfa/htpolynet-tests/testset/3-bisgma-styrene-thermoset/proj-0
    INFO> *************************** Ambertools: ***************************
    INFO>   antechamber (ver. 22.0) at antechamber                                        
    INFO>         tleap (ver. 22.0) at tleap                                              
    INFO>      parmchk2 (ver. 22.0) at parmchk2                                           
    INFO> Configuration: GMASTY-lo.yaml

Then a description of all the molecules that need to be created:

.. code-block:: console

    INFO> *********** Templates in proj-0/molecules/parameterized ***********
    INFO> 32 molecules detected in GMASTY-lo.yaml
    INFO>                       explicit: 11   
    INFO>     implied by stereochemistry: 21   
    INFO>            implied by symmetry: 0    

First are the molecules that are explicit in the configuration file:

.. code-block:: console

    INFO> AmberTools> generating GAFF parameters from BPA.mol2
    INFO> BPA: 228.28 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE.mol2
    INFO> HIE: 146.18 g/mol
    INFO> AmberTools> generating GAFF parameters from STY.mol2
    INFO> STY: 106.16 g/mol
    INFO> AmberTools> generating GAFF parameters from GM1.mol2
    INFO> GM1: 372.44 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~STY.mol2
    INFO> STY~C1-C2~STY: 210.30 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~HIE.mol2
    INFO> HIE~C1-C2~HIE: 290.35 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~HIE.mol2
    INFO> STY~C1-C2~HIE: 250.33 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~STY.mol2
    INFO> HIE~C1-C2~STY: 250.33 g/mol
    INFO> AmberTools> generating GAFF parameters from STYCC.mol2
    INFO> STYCC: 104.14 g/mol
    INFO> AmberTools> generating GAFF parameters from HIECC.mol2
    INFO> HIECC: 144.17 g/mol
    INFO> AmberTools> generating GAFF parameters from GMA.mol2
    INFO> GMA: 516.61 g/mol

Next are the diastereomers that are implied by stereocenters.  Note that these are just built, not parameterized:

.. code-block:: console

    INFO> Built GM1S-1 using topology of GM1; copying GM1.top to GM1S-1.top
    INFO> GM1S-1: 372.44 g/mol
    INFO> Built GM1S-2 using topology of GM1; copying GM1.top to GM1S-2.top
    INFO> GM1S-2: 372.44 g/mol
    INFO> Built GM1S-3 using topology of GM1; copying GM1.top to GM1S-3.top
    INFO> GM1S-3: 372.44 g/mol
    INFO> Built GMAS-1 using topology of GMA; copying GMA.top to GMAS-1.top
    INFO> GMAS-1: 516.61 g/mol
    INFO> Built GMAS-2 using topology of GMA; copying GMA.top to GMAS-2.top
    INFO> GMAS-2: 516.61 g/mol
    INFO> Built GMAS-3 using topology of GMA; copying GMA.top to GMAS-3.top
    INFO> GMAS-3: 516.61 g/mol
    INFO> Built GMAS-4 using topology of GMA; copying GMA.top to GMAS-4.top
    INFO> GMAS-4: 516.61 g/mol
    INFO> Built GMAS-5 using topology of GMA; copying GMA.top to GMAS-5.top
    INFO> GMAS-5: 516.61 g/mol
    INFO> Built GMAS-6 using topology of GMA; copying GMA.top to GMAS-6.top
    INFO> GMAS-6: 516.61 g/mol
    INFO> Built GMAS-7 using topology of GMA; copying GMA.top to GMAS-7.top
    INFO> GMAS-7: 516.61 g/mol
    INFO> Built GMAS-8 using topology of GMA; copying GMA.top to GMAS-8.top
    INFO> GMAS-8: 516.61 g/mol
    INFO> Built GMAS-9 using topology of GMA; copying GMA.top to GMAS-9.top
    INFO> GMAS-9: 516.61 g/mol
    INFO> Built GMAS-10 using topology of GMA; copying GMA.top to GMAS-10.top
    INFO> GMAS-10: 516.61 g/mol
    INFO> Built GMAS-11 using topology of GMA; copying GMA.top to GMAS-11.top
    INFO> GMAS-11: 516.61 g/mol
    INFO> Built GMAS-12 using topology of GMA; copying GMA.top to GMAS-12.top
    INFO> GMAS-12: 516.61 g/mol
    INFO> Built GMAS-13 using topology of GMA; copying GMA.top to GMAS-13.top
    INFO> GMAS-13: 516.61 g/mol
    INFO> Built GMAS-14 using topology of GMA; copying GMA.top to GMAS-14.top
    INFO> GMAS-14: 516.61 g/mol
    INFO> Built GMAS-15 using topology of GMA; copying GMA.top to GMAS-15.top
    INFO> GMAS-15: 516.61 g/mol

Notice that this generates three of the four diastereomers of the intermediate GM1, and 15 of the 16 diastereomers of GMA; the missing ones were already generated!  This set of 16 GMA diastereomers is used to add GMA molecules to the initial liquid.

Finally, the 32 molecules implied by chain-expansion of the cure reactions:

.. code-block:: console

    INFO> 32 molecules implied by chaining
    INFO> AmberTools> generating GAFF parameters from HIE~C1=C2~STY~C1-C2~STY.mol2
    INFO> HIE~C1=C2~STY~C1-C2~STY: 354.47 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1=C2~HIE~C1-C2~HIE.mol2
    INFO> HIE~C1=C2~HIE~C1-C2~HIE: 434.51 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1=C2~STY~C1-C2~HIE.mol2
    INFO> HIE~C1=C2~STY~C1-C2~HIE: 394.49 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1=C2~HIE~C1-C2~STY.mol2
    INFO> HIE~C1=C2~HIE~C1-C2~STY: 394.49 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1=C2~STY~C1-C2~STY.mol2
    INFO> STY~C1=C2~STY~C1-C2~STY: 314.45 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1=C2~HIE~C1-C2~HIE.mol2
    INFO> STY~C1=C2~HIE~C1-C2~HIE: 394.49 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1=C2~STY~C1-C2~HIE.mol2
    INFO> STY~C1=C2~STY~C1-C2~HIE: 354.47 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1=C2~HIE~C1-C2~STY.mol2
    INFO> STY~C1=C2~HIE~C1-C2~STY: 354.47 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~STY~C1=C2~HIE.mol2
    INFO> STY~C1-C2~STY~C1=C2~HIE: 354.47 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~HIE~C1=C2~HIE.mol2
    INFO> HIE~C1-C2~HIE~C1=C2~HIE: 434.51 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~HIE~C1=C2~HIE.mol2
    INFO> STY~C1-C2~HIE~C1=C2~HIE: 394.49 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~STY~C1=C2~HIE.mol2
    INFO> HIE~C1-C2~STY~C1=C2~HIE: 394.49 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~STY~C1=C2~STY.mol2
    INFO> STY~C1-C2~STY~C1=C2~STY: 314.45 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~HIE~C1=C2~STY.mol2
    INFO> HIE~C1-C2~HIE~C1=C2~STY: 394.49 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~HIE~C1=C2~STY.mol2
    INFO> STY~C1-C2~HIE~C1=C2~STY: 354.47 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~STY~C1=C2~STY.mol2
    INFO> HIE~C1-C2~STY~C1=C2~STY: 354.47 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~STY~C1=C2~STY~C1-C2~STY.mol2
    INFO> STY~C1-C2~STY~C1=C2~STY~C1-C2~STY: 418.59 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~STY~C1=C2~HIE~C1-C2~HIE.mol2
    INFO> STY~C1-C2~STY~C1=C2~HIE~C1-C2~HIE: 498.64 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~STY~C1=C2~STY~C1-C2~HIE.mol2
    INFO> STY~C1-C2~STY~C1=C2~STY~C1-C2~HIE: 458.61 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~STY~C1=C2~HIE~C1-C2~STY.mol2
    INFO> STY~C1-C2~STY~C1=C2~HIE~C1-C2~STY: 458.61 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~HIE~C1=C2~STY~C1-C2~STY.mol2
    INFO> HIE~C1-C2~HIE~C1=C2~STY~C1-C2~STY: 498.64 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~HIE~C1=C2~HIE~C1-C2~HIE.mol2
    INFO> HIE~C1-C2~HIE~C1=C2~HIE~C1-C2~HIE: 578.68 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~HIE~C1=C2~STY~C1-C2~HIE.mol2
    INFO> HIE~C1-C2~HIE~C1=C2~STY~C1-C2~HIE: 538.66 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~HIE~C1=C2~HIE~C1-C2~STY.mol2
    INFO> HIE~C1-C2~HIE~C1=C2~HIE~C1-C2~STY: 538.66 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~HIE~C1=C2~STY~C1-C2~STY.mol2
    INFO> STY~C1-C2~HIE~C1=C2~STY~C1-C2~STY: 458.61 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~HIE~C1=C2~HIE~C1-C2~HIE.mol2
    INFO> STY~C1-C2~HIE~C1=C2~HIE~C1-C2~HIE: 538.66 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~HIE~C1=C2~STY~C1-C2~HIE.mol2
    INFO> STY~C1-C2~HIE~C1=C2~STY~C1-C2~HIE: 498.64 g/mol
    INFO> AmberTools> generating GAFF parameters from STY~C1-C2~HIE~C1=C2~HIE~C1-C2~STY.mol2
    INFO> STY~C1-C2~HIE~C1=C2~HIE~C1-C2~STY: 498.64 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~STY~C1=C2~STY~C1-C2~STY.mol2
    INFO> HIE~C1-C2~STY~C1=C2~STY~C1-C2~STY: 458.61 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~STY~C1=C2~HIE~C1-C2~HIE.mol2
    INFO> HIE~C1-C2~STY~C1=C2~HIE~C1-C2~HIE: 538.66 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~STY~C1=C2~STY~C1-C2~HIE.mol2
    INFO> HIE~C1-C2~STY~C1=C2~STY~C1-C2~HIE: 498.64 g/mol
    INFO> AmberTools> generating GAFF parameters from HIE~C1-C2~STY~C1=C2~HIE~C1-C2~STY.mol2
    INFO> HIE~C1-C2~STY~C1=C2~HIE~C1-C2~STY: 498.64 g/mol
    INFO> Generated 64 molecule templates
    INFO> Initial composition is STY 150, GMA 75
    INFO> 100% conversion is 300 bonds


Now the creation and densification of the initial liquid:

.. code-block:: console

    INFO> ************** Initialization in proj-0/systems/init **************
    INFO> Topology "init.top" in proj-0/systems/init
    INFO> Initial density: 100.0 kg/m^3
    INFO> Total mass: 9.078e-23 kg
    INFO> Box aspect ratio: 1.0 x 1.0 x 1.0
    INFO> Initial box side lengths: 9.683 nm x 9.683 nm x 9.683 nm
    INFO> Coordinates "init.gro" in proj-0/systems/init
    INFO> Extended attributes "init.grx" in proj-0/systems/init
    INFO> ********** Densification in proj-0/systems/densification **********
    INFO> Running Gromacs: minimization
    INFO> Running Gromacs: nvt ensemble;  10.00 ps,  300.00 K
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,  10.00 bar
    INFO> Current box side lengths: 9.256 nm x 9.256 nm x 9.256 nm
    INFO> Density                      114.46
    INFO> Running-average-Density      106.81
    INFO> Rolling-10-average-Density   113.66
    INFO> Repeat 1 out of 8
    INFO> Current box side lengths: 8.773 nm x 8.773 nm x 8.773 nm
    INFO> Density                      134.47
    INFO> Running-average-Density      123.99
    INFO> Rolling-10-average-Density   133.40
    INFO> Repeat 2 out of 8
    INFO> Current box side lengths: 8.242 nm x 8.242 nm x 8.242 nm
    INFO> Density                      162.15
    INFO> Running-average-Density      147.43
    INFO> Rolling-10-average-Density   160.62
    INFO> Repeat 3 out of 8
    INFO> Current box side lengths: 7.573 nm x 7.573 nm x 7.573 nm
    INFO> Density                      209.06
    INFO> Running-average-Density      183.17
    INFO> Rolling-10-average-Density   206.17
    INFO> Repeat 4 out of 8
    INFO> Current box side lengths: 6.589 nm x 6.589 nm x 6.589 nm
    INFO> Density                      317.36
    INFO> Running-average-Density      254.44
    INFO> Rolling-10-average-Density   310.28
    INFO> Repeat 5 out of 8
    INFO> Current box side lengths: 4.569 nm x 4.569 nm x 4.569 nm
    INFO> Density                      952.04
    INFO> Running-average-Density      547.44
    INFO> Rolling-10-average-Density   921.37
    INFO> Repeat 6 out of 8
    INFO> Current box side lengths: 4.482 nm x 4.482 nm x 4.482 nm
    INFO> Density                      1008.60
    INFO> Running-average-Density      1001.76
    INFO> Rolling-10-average-Density   1006.40
    INFO> Repeat 7 out of 8
    INFO> Current box side lengths: 4.475 nm x 4.475 nm x 4.475 nm
    INFO> Density                      1013.33
    INFO> Running-average-Density      1013.13
    INFO> Rolling-10-average-Density   1017.67
    INFO> Repeat 8 out of 8
    INFO> Current box side lengths: 4.469 nm x 4.469 nm x 4.469 nm
    INFO> Density                      1016.87
    INFO> Running-average-Density      1016.31
    INFO> Rolling-10-average-Density   1018.33
    INFO> Densified coordinates in proj-0/systems/densification/densified-repeat-8-npt.gro

Now the precure equilibration and annealing:

.. code-block:: console

    INFO> **************** Precure in proj-0/systems/precure ****************
    INFO> Running Gromacs: npt ensemble; 300.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.465 nm x 4.465 nm x 4.465 nm
    INFO> Density                      1019.94
    INFO> Running-average-Density      1018.01
    INFO> Rolling-10-average-Density   1019.32
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 200.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.463 nm x 4.463 nm x 4.463 nm
    INFO> Density                      1021.26
    INFO> Running-average-Density      1016.82
    INFO> Rolling-10-average-Density   1020.27


Note that we achieve a pretty good initial density for this liquid of about 1.02 g/cc.  Now the CURE algorithm begins, seeking to form 150 out of the total possible 300 bonds:

.. code-block:: console

    INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> Attempting to form 150 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 1 will generate 77 new bonds
    INFO> Cure_relax initiated on 77 distances (max 0.483 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.426                  0.651
    INFO>          2              0.385                  0.611
    INFO>          3              0.330                  0.551
    INFO>          4              0.278                  0.507
    INFO>          5              0.215                  0.464
    INFO>          6              0.166                  0.419
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.424 nm x 4.424 nm x 4.424 nm
    INFO> Density                      1045.51
    INFO> Running-average-Density      1012.53
    INFO> Rolling-10-average-Density   1046.27
    INFO> Iteration 1 current conversion 0.257 or 77 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 64 new bonds
    INFO> Cure_relax initiated on 64 distances (max 0.492 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.431                  0.665
    INFO>          2              0.386                  0.617
    INFO>          3              0.328                  0.565
    INFO>          4              0.272                  0.507
    INFO>          5              0.214                  0.469
    INFO>          6              0.169                  0.430
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.391 nm x 4.391 nm x 4.391 nm
    INFO> Density                      1066.64
    INFO> Running-average-Density      1038.46
    INFO> Rolling-10-average-Density   1065.77
    INFO> Iteration 2 current conversion 0.470 or 141 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 3 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 3 will generate 10 new bonds
    INFO> Cure_relax initiated on 10 distances (max 0.410 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.365                  0.596
    INFO>          2              0.315                  0.544
    INFO>          3              0.261                  0.513
    INFO>          4              0.221                  0.452
    INFO>          5              0.166                  0.410
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.390 nm x 4.390 nm x 4.390 nm
    INFO> Density                      1066.85
    INFO> Running-average-Density      1044.10
    INFO> Rolling-10-average-Density   1070.23
    INFO> Iteration 3 current conversion 0.503 or 151 bonds

We see here that only three CURE iterations were required.  Now comes capping of the 55 activated doule-bonds that neither attacked nor were attacked:

.. code-block:: console

    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 55 new bonds
    INFO> Cap_relax initiated on 55 distances (max 0.155 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.145                  0.382
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.367 nm x 4.367 nm x 4.367 nm
    INFO> Density                      1082.01
    INFO> Running-average-Density      1074.59
    INFO> Rolling-10-average-Density   1082.57
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********

And finally, postcure equilibration and annealing:

.. code-block:: console

    INFO> *************** Postcure in proj-0/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.357 nm x 4.357 nm x 4.357 nm
    INFO> Density                      1089.52
    INFO> Running-average-Density      1087.44
    INFO> Rolling-10-average-Density   1087.69
    INFO> *********** Final data to proj-0/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

So we see our 50% cured system reached about 1.09 g/cc.  I am not claiming any of these equilibration runs are long enough or the system is big enough, but this example just serves to illustrate how ``HTPolyNet`` works.

If we consider the ``hi``-cure system, we see that it took 19 iterations to reach 95% cure.  Below is the end of the console output for that run:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 19 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Iteration 19 will generate 1 new bond
    INFO> Cure_relax initiated on 1 distance (max 0.554 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.487                  0.692
    INFO>          2              0.420                  0.607
    INFO>          3              0.355                  0.569
    INFO>          4              0.284                  0.503
    INFO>          5              0.230                  0.446
    INFO>          6              0.159                  0.408
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.398 nm x 4.398 nm x 4.398 nm
    INFO> Density                      1055.67
    INFO> Running-average-Density      1039.04
    INFO> Rolling-10-average-Density   1052.93
    INFO> Iteration 19 current conversion 0.950 or 285 bonds
    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 0 new bonds
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********
    INFO> *************** Postcure in proj-1/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.366 nm x 4.366 nm x 4.366 nm
    INFO> Density                      1079.18
    INFO> Running-average-Density      1077.62
    INFO> Rolling-10-average-Density   1080.11
    INFO> *********** Final data to proj-1/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

Notice that no capping reactions happen since all activated double-bonds are either attacked or attackers or both.

Let's now take a look at some selected :ref:`results <ve_results>`.