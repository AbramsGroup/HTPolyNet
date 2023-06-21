.. _ve_run:

Running the Build
=================

Now, in our working directory ``my_vinyl_ester/3-bisgma-styrene-thermoset/``, if you haven't already, go ahead and launch via ``./run.sh``.

.. code-block:: console

    $ ./run.sh

As the two ``htpolynet run`` commands indicate, standard output is being redirected from the console to the file ``console.log``:

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
    INFO> User library is /home/cfa/htpolynet-tutorials/v1.0.6/3-bisgma-styrene-thermoset/lib
    INFO> New project in /home/cfa/htpolynet-tutorials/v1.0.6/3-bisgma-styrene-thermoset/proj-0
    INFO> *************************** Ambertools: ***************************
    INFO> ********************  antechamber (ver. 22.0) *********************
    INFO> ********************        tleap (ver. 22.0) *********************
    INFO> ********************     parmchk2 (ver. 22.0) *********************
    INFO> Configuration: GMASTY.yaml

Then a description of all the molecules that need to be created:

.. code-block:: console

    INFO> *********** Templates in proj-0/molecules/parameterized ***********
    INFO> 32 molecules detected in GMASTY.yaml
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
    INFO> Current box side lengths: 9.240 nm x 9.240 nm x 9.240 nm
    INFO> Density                      115.07
    INFO> Running-average-Density      107.18
    INFO> Rolling-10-average-Density   114.32
    INFO> Repeat 1 out of 8
    INFO> Current box side lengths: 8.761 nm x 8.761 nm x 8.761 nm
    INFO> Density                      135.00
    INFO> Running-average-Density      124.69
    INFO> Rolling-10-average-Density   133.93
    INFO> Repeat 2 out of 8
    INFO> Current box side lengths: 8.215 nm x 8.215 nm x 8.215 nm
    INFO> Density                      163.72
    INFO> Running-average-Density      148.41
    INFO> Rolling-10-average-Density   162.23
    INFO> Repeat 3 out of 8
    INFO> Current box side lengths: 7.565 nm x 7.565 nm x 7.565 nm
    INFO> Density                      209.65
    INFO> Running-average-Density      184.25
    INFO> Rolling-10-average-Density   206.68
    INFO> Repeat 4 out of 8
    INFO> Current box side lengths: 6.544 nm x 6.544 nm x 6.544 nm
    INFO> Density                      323.89
    INFO> Running-average-Density      257.92
    INFO> Rolling-10-average-Density   316.64
    INFO> Repeat 5 out of 8
    INFO> Current box side lengths: 4.533 nm x 4.533 nm x 4.533 nm
    INFO> Density                      974.55
    INFO> Running-average-Density      572.15
    INFO> Rolling-10-average-Density   949.96
    INFO> Repeat 6 out of 8
    INFO> Current box side lengths: 4.475 nm x 4.475 nm x 4.475 nm
    INFO> Density                      1012.95
    INFO> Running-average-Density      1007.48
    INFO> Rolling-10-average-Density   1015.49
    INFO> Repeat 7 out of 8
    INFO> Current box side lengths: 4.473 nm x 4.473 nm x 4.473 nm
    INFO> Density                      1014.68
    INFO> Running-average-Density      1014.93
    INFO> Rolling-10-average-Density   1016.80
    INFO> Repeat 8 out of 8
    INFO> Current box side lengths: 4.463 nm x 4.463 nm x 4.463 nm
    INFO> Density                      1021.02
    INFO> Running-average-Density      1019.14
    INFO> Rolling-10-average-Density   1022.00
    INFO> Densified coordinates in proj-0/systems/densification/densified-repeat-8-npt.gro

Note that we achieve a pretty good initial density for this liquid of about 1.02 g/cc.  Now the CURE algorithm begins, seeking to form 285 out of the total possible 300 bonds.  The console output shows that 161 bonds form in just the first two iterations:

.. code-block:: console

    INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> Attempting to form 285 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 1 will generate 90 new bonds
    INFO> Step "cure_relax" initiated on 90 distances (max 0.489 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.440                  0.672
    INFO>          2              0.387                  0.602
    INFO>          3              0.328                  0.561
    INFO>          4              0.267                  0.508
    INFO>          5              0.216                  0.458
    INFO>          6              0.168                  0.418
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.417 nm x 4.417 nm x 4.417 nm
    INFO> Density                      1049.89
    INFO> Running-average-Density      1021.21
    INFO> Rolling-10-average-Density   1049.73
    INFO> Iteration 1 current conversion 0.300 or 90 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 71 new bonds
    INFO> Step "cure_relax" initiated on 71 distances (max 0.499 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.434                  0.656
    INFO>          2              0.389                  0.625
    INFO>          3              0.334                  0.574
    INFO>          4              0.271                  0.516
    INFO>          5              0.215                  0.459
    INFO>          6              0.169                  0.426
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.384 nm x 4.384 nm x 4.384 nm
    INFO> Density                      1071.08
    INFO> Running-average-Density      1038.69
    INFO> Rolling-10-average-Density   1067.55
    INFO> Iteration 2 current conversion 0.537 or 161 bonds

For this build, a total of 22 CURE iterations were required to reach 95% conversion.  Here is the console output for the last two iterations, and we see they each generated only one bond, and both required pre-bond dragging before bond formation and relaxation:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 21 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Radius increased to 1.0 nm
    INFO> Iteration 21 will generate 1 new bond
    INFO> Step "cure_drag" initiated on 1 distance (max 0.961 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.942
    INFO>          2              0.854
    INFO>          3              0.786
    INFO>          4              0.743
    INFO>          5              0.689
    INFO>          6              0.631
    INFO>          7              0.576
    INFO>          8              0.523
    INFO>          9              0.462
    INFO>         10              0.418
    INFO>         11              0.359
    INFO>         12              0.310
    INFO> Step "cure_relax" initiated on 1 distance (max 0.310 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.262                  0.499
    INFO>          2              0.218                  0.461
    INFO>          3              0.161                  0.405
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.393 nm x 4.393 nm x 4.393 nm
    INFO> Density                      1059.61
    INFO> Running-average-Density      1049.49
    INFO> Rolling-10-average-Density   1062.18
    INFO> Iteration 21 current conversion 0.947 or 284 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 22 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Radius increased to 1.0 nm
    INFO> Iteration 22 will generate 1 new bond
    INFO> Step "cure_drag" initiated on 1 distance (max 0.977 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.926
    INFO>          2              0.860
    INFO>          3              0.807
    INFO>          4              0.747
    INFO>          5              0.691
    INFO>          6              0.631
    INFO>          7              0.580
    INFO>          8              0.533
    INFO>          9              0.475
    INFO>         10              0.424
    INFO>         11              0.359
    INFO>         12              0.308
    INFO> Step "cure_relax" initiated on 1 distance (max 0.308 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.274                  0.508
    INFO>          2              0.208                  0.438
    INFO>          3              0.161                  0.415
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.389 nm x 4.389 nm x 4.389 nm
    INFO> Density                      1062.66
    INFO> Running-average-Density      1047.72
    INFO> Rolling-10-average-Density   1060.33
    INFO> Iteration 22 current conversion 0.950 or 285 bonds

And since there are apparently no unreactd double bonds, there are no capping reactions necessary:

.. code-block:: console

    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 0 new bonds
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********

And finally, postcure equilibration and annealing:

.. code-block:: console

    INFO> *************** Postcure in proj-0/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.362 nm x 4.362 nm x 4.362 nm
    INFO> Density                      1082.56
    INFO> Running-average-Density      1080.87
    INFO> Rolling-10-average-Density   1083.25
    INFO> *********** Final data to proj-0/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

So we see our 95% cured system reached about 1.08 g/cc.  I am not claiming any of these equilibration runs are long enough or the system is big enough, but this example just serves to illustrate how ``HTPolyNet`` works.

Let's now take a look at some selected :ref:`results <ve_results>`.