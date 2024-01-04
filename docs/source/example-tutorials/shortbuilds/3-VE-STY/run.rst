.. _ve_run:

Running the Build
=================

Now, in our working directory ``my_vinyl_ester/3-bisgma-styrene-thermoset/``, if you haven't already, go ahead and launch via ``./run.sh``.

.. code-block:: console

    $ ./run.sh

As the two ``htpolynet run`` commands indicate, standard output is being redirected from the console to the file ``console.log``:

.. code-block:: console

    INFO>                                                                    
    INFO>     HTPolyNet 1.0.8                                                
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
    INFO>     Please cite the HTPolyNet paper:                               
    INFO>                                                                    
    INFO>     Ming Huang and Cameron F. Abrams, HTPolyNet: A general         
    INFO>     system generator for all-atom molecular simulations of         
    INFO>     amorphous crosslinked polymers, SoftwareX, vol. 21,            
    INFO>     pp. 101303, 2023 (doi:10.1016/j.softx.2022.101303) 
    INFO>                                                                    
    INFO> ******************** HTPolyNet runtime begins *********************
    INFO> User library is /home/cfa/htpolynet-tutorials/1.0.8/3-bisgma-styrene-thermoset/lib
    INFO> New project in /home/cfa/htpolynet-tutorials/1.0.8/3-bisgma-styrene-thermoset/proj-0
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
    INFO> Current box side lengths: 9.215 nm x 9.215 nm x 9.215 nm
    INFO> Density                      116.01
    INFO> Running-average-Density      107.13
    INFO> Rolling-10-average-Density   114.85
    INFO> Repeat 1 out of 8
    INFO> Current box side lengths: 8.634 nm x 8.634 nm x 8.634 nm
    INFO> Density                      141.04
    INFO> Running-average-Density      127.82
    INFO> Rolling-10-average-Density   139.92
    INFO> Repeat 2 out of 8
    INFO> Current box side lengths: 7.969 nm x 7.969 nm x 7.969 nm
    INFO> Density                      179.39
    INFO> Running-average-Density      157.97
    INFO> Rolling-10-average-Density   177.29
    INFO> Repeat 3 out of 8
    INFO> Current box side lengths: 7.283 nm x 7.283 nm x 7.283 nm
    INFO> Density                      235.01
    INFO> Running-average-Density      204.78
    INFO> Rolling-10-average-Density   231.27
    INFO> Repeat 4 out of 8
    INFO> Current box side lengths: 6.375 nm x 6.375 nm x 6.375 nm
    INFO> Density                      350.34
    INFO> Running-average-Density      279.94
    INFO> Rolling-10-average-Density   341.25
    INFO> Repeat 5 out of 8
    INFO> Current box side lengths: 4.669 nm x 4.669 nm x 4.669 nm
    INFO> Density                      891.80
    INFO> Running-average-Density      548.19
    INFO> Rolling-10-average-Density   867.14
    INFO> Repeat 6 out of 8
    INFO> Current box side lengths: 4.462 nm x 4.462 nm x 4.462 nm
    INFO> Density                      1021.60
    INFO> Running-average-Density       993.90
    INFO> Rolling-10-average-Density   1015.59
    INFO> Repeat 7 out of 8
    INFO> Current box side lengths: 4.471 nm x 4.471 nm x 4.471 nm
    INFO> Density                      1015.42
    INFO> Running-average-Density      1013.50
    INFO> Rolling-10-average-Density   1012.53
    INFO> Repeat 8 out of 8
    INFO> Current box side lengths: 4.465 nm x 4.465 nm x 4.465 nm
    INFO> Density                      1020.05
    INFO> Running-average-Density      1016.22
    INFO> Rolling-10-average-Density   1016.10
    INFO> Densified coordinates in proj-0/systems/densification/densified-repeat-8-npt.gro

Note that we achieve a pretty good initial density for this liquid of about 1.02 g/cc.  Now the CURE algorithm begins, seeking to form 285 out of the total possible 300 bonds.  The console output shows that 159 bonds form in just the first two iterations:

.. code-block:: console

    INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> Attempting to form 285 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 1 will generate 81 new bonds
    INFO> Step "cure_relax" initiated on 81 distances (max 0.493 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.454                  0.660
    INFO>          2              0.386                  0.607
    INFO>          3              0.328                  0.556
    INFO>          4              0.275                  0.504
    INFO>          5              0.215                  0.461
    INFO>          6              0.168                  0.417
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.424 nm x 4.424 nm x 4.424 nm
    INFO> Density                      1045.35
    INFO> Running-average-Density      1013.18
    INFO> Rolling-10-average-Density   1044.08
    INFO> Iteration 1 current conversion 0.270 or 81 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 78 new bonds
    INFO> Step "cure_relax" initiated on 78 distances (max 0.499 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.446                  0.668
    INFO>          2              0.394                  0.635
    INFO>          3              0.326                  0.565
    INFO>          4              0.272                  0.517
    INFO>          5              0.218                  0.465
    INFO>          6              0.167                  0.423
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.393 nm x 4.393 nm x 4.393 nm
    INFO> Density                      1064.22
    INFO> Running-average-Density      1038.23
    INFO> Rolling-10-average-Density   1061.93
    INFO> Iteration 2 current conversion 0.530 or 159 bonds

For this build, a total of 16 CURE iterations were required to reach 95% conversion.  Here is the console output for the last two iterations, and we see they generated one and two bonds, respectively, and the last one required pre-bond dragging before bond formation and relaxation:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 15 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Iteration 15 will generate 1 new bond
    INFO> Step "cure_relax" initiated on 1 distance (max 0.543 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.475                  0.674
    INFO>          2              0.417                  0.634
    INFO>          3              0.349                  0.559
    INFO>          4              0.286                  0.518
    INFO>          5              0.229                  0.479
    INFO>          6              0.158                  0.428
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.393 nm x 4.393 nm x 4.393 nm
    INFO> Density                      1059.35
    INFO> Running-average-Density      1045.04
    INFO> Rolling-10-average-Density   1060.65
    INFO> Iteration 15 current conversion 0.943 or 283 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 16 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Iteration 16 will generate 2 new bonds
    INFO> Step "cure_drag" initiated on 2 distances (max 0.691 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.645
    INFO>          2              0.604
    INFO>          3              0.555
    INFO>          4              0.497
    INFO>          5              0.445
    INFO>          6              0.398
    INFO>          7              0.356
    INFO>          8              0.306
    INFO> Step "cure_relax" initiated on 2 distances (max 0.306 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.259                  0.501
    INFO>          2              0.210                  0.451
    INFO>          3              0.165                  0.421
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 4.396 nm x 4.396 nm x 4.396 nm
    INFO> Density                      1057.45
    INFO> Running-average-Density      1049.88
    INFO> Rolling-10-average-Density   1056.57
    INFO> Iteration 16 current conversion 0.950 or 285 bonds

And since there are apparently no unreacted double bonds, there are no capping reactions necessary:

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
    INFO> Current box side lengths: 4.358 nm x 4.358 nm x 4.358 nm
    INFO> Density                      1085.65
    INFO> Running-average-Density      1083.04
    INFO> Rolling-10-average-Density   1082.58
    INFO> *********** Final data to proj-0/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

So we see our 95% cured system reached about 1.08 g/cc.  I am not claiming any of these equilibration runs are long enough or the system is big enough, but this example just serves to illustrate how ``HTPolyNet`` works.

Let's now take a look at some selected :ref:`results <ve_results>`.