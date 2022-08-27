.. _dgeba_run:

Running the Build
=================

Now, in our working directory ``my_dgeba_pacm/4-pacm-dgeba-epoxy-thermoset/``, if you haven't already, go ahead and launch via ``./run.sh``.

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
    INFO> User library is /home/cfa/htpolynet-tests/testset/4-pacm-dgeba-epoxy-thermoset/lib
    INFO> New project in /home/cfa/htpolynet-tests/testset/4-pacm-dgeba-epoxy-thermoset/proj-0
    INFO> *************************** Ambertools: ***************************
    INFO>   antechamber (ver. 22.0) at antechamber                                        
    INFO>         tleap (ver. 22.0) at tleap                                              
    INFO>      parmchk2 (ver. 22.0) at parmchk2                                           
    INFO> Configuration: DGE-PAC-lo.yaml

Then a description of all the molecules that need to be created:

.. code-block:: console

    INFO> *********** Templates in proj-0/molecules/parameterized ***********
    INFO> 22 molecules detected in DGE-PAC-lo.yaml
    INFO>                       explicit: 5    
    INFO>     implied by stereochemistry: 6    
    INFO>            implied by symmetry: 11   
    INFO> AmberTools> generating GAFF parameters from PAC.mol2
    INFO> PAC: 210.36 g/mol
    INFO> AmberTools> generating GAFF parameters from DGE.mol2
    INFO> DGE: 344.43 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N1-C1~DGE.mol2
    INFO> PAC~N1-C1~DGE: 552.78 g/mol
    INFO> AmberTools> generating GAFF parameters from DGEC.mol2
    INFO> DGEC: 342.42 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N1-C1~DGE-C1~DGE.mol2
    INFO> PAC~N1-C1~DGE-C1~DGE: 895.19 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N1-C2~DGE.mol2
    INFO> PAC~N1-C2~DGE: 552.78 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N1-C2~DGE~N1-C1~DGE.mol2
    INFO> PAC~N1-C2~DGE~N1-C1~DGE: 895.19 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N2-C1~DGE.mol2
    INFO> PAC~N2-C1~DGE: 552.78 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N2-C1~DGE~N2-C1~DGE.mol2
    INFO> PAC~N2-C1~DGE~N2-C1~DGE: 895.19 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N2-C2~DGE.mol2
    INFO> PAC~N2-C2~DGE: 552.78 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N2-C2~DGE~N2-C1~DGE.mol2
    INFO> PAC~N2-C2~DGE~N2-C1~DGE: 895.19 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N1-C1~DGE~N1-C2~DGE.mol2
    INFO> PAC~N1-C1~DGE~N1-C2~DGE: 895.19 g/mol
    INFO> AmberTools> generating GAFF parameters from DGEC-1.mol2
    INFO> DGEC-1: 342.42 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N1-C2~DGE~N1-C2~DGE.mol2
    INFO> PAC~N1-C2~DGE~N1-C2~DGE: 895.19 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N2-C1~DGE~N2-C2~DGE.mol2
    INFO> PAC~N2-C1~DGE~N2-C2~DGE: 895.19 g/mol
    INFO> AmberTools> generating GAFF parameters from PAC~N2-C2~DGE~N2-C2~DGE.mol2
    INFO> PAC~N2-C2~DGE~N2-C2~DGE: 895.19 g/mol
    INFO> Generated 22 molecule templates
    INFO> Initial composition is DGE 200, PAC 100
    INFO> 100% conversion is 400 bonds

Notice that the symmetry-equivalent atoms in both monomers result in a total of 16 unique cure reactions forming secondary and tertiary amines, as well as two capping reactions for DGEBA.  Now the initial construction of the system and its densification are next:

.. code-block:: console

    INFO> ************** Initialization in proj-0/systems/init **************
    INFO> Topology "init.top" in proj-0/systems/init
    INFO> Initial density: 300.0 kg/m^3
    INFO> Total mass: 1.493e-22 kg
    INFO> Box aspect ratio: 1.0 x 1.0 x 1.0
    INFO> Initial box side lengths: 7.925 nm x 7.925 nm x 7.925 nm
    INFO> Coordinates "init.gro" in proj-0/systems/init
    INFO> Extended attributes "init.grx" in proj-0/systems/init
    INFO> ********** Densification in proj-0/systems/densification **********
    INFO> Running Gromacs: minimization
    INFO> Running Gromacs: nvt ensemble;  10.00 ps,  300.00 K
    INFO> Running Gromacs: npt ensemble; 300.00 ps,  300.00 K,  10.00 bar
    INFO> Current box side lengths: 5.228 nm x 5.228 nm x 5.228 nm
    INFO> Density                      1045.16
    INFO> Running-average-Density       893.88
    INFO> Rolling-10-average-Density   1045.73
    INFO> Densified coordinates in proj-0/systems/densification/densified-npt.gro
    INFO> **************** Precure in proj-0/systems/precure ****************
    INFO> Running Gromacs: npt ensemble; 200.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.227 nm x 5.227 nm x 5.227 nm
    INFO> Density                      1045.38
    INFO> Running-average-Density      1048.43
    INFO> Rolling-10-average-Density   1047.02
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.236 nm x 5.236 nm x 5.236 nm
    INFO> Density                      1039.97
    INFO> Running-average-Density      1041.91
    INFO> Rolling-10-average-Density   1040.19

Note that we achieve a pretty good initial density for this liquid of about 1.04 g/cc.  Now the CURE algorithm begins, seeking to form 200 out of the total possible 400 bonds:

.. code-block:: console

    INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> Attempting to form 200 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 1 will generate 71 new bonds
    INFO> Cure_relax initiated on 71 distances (max 0.485 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.439                  0.665
    INFO>          2              0.373                  0.591
    INFO>          3              0.311                  0.541
    INFO>          4              0.270                  0.494
    INFO>          5              0.207                  0.455
    INFO>          6              0.156                  0.404
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.227 nm x 5.227 nm x 5.227 nm
    INFO> Density                      1043.88
    INFO> Running-average-Density      1020.22
    INFO> Rolling-10-average-Density   1043.70
    INFO> Iteration 1 current conversion 0.177 or 71 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 55 new bonds
    INFO> Cure_relax initiated on 55 distances (max 0.485 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.441                  0.670
    INFO>          2              0.370                  0.599
    INFO>          3              0.313                  0.540
    INFO>          4              0.264                  0.500
    INFO>          5              0.204                  0.448
    INFO>          6              0.157                  0.415
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.214 nm x 5.214 nm x 5.214 nm
    INFO> Density                      1050.68
    INFO> Running-average-Density      1028.62
    INFO> Rolling-10-average-Density   1050.15
    INFO> Iteration 2 current conversion 0.315 or 126 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 3 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 3 will generate 36 new bonds
    INFO> Cure_relax initiated on 36 distances (max 0.464 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.397                  0.632
    INFO>          2              0.349                  0.571
    INFO>          3              0.283                  0.529
    INFO>          4              0.214                  0.457
    INFO>          5              0.156                  0.413
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.204 nm x 5.204 nm x 5.204 nm
    INFO> Density                      1055.60
    INFO> Running-average-Density      1033.13
    INFO> Rolling-10-average-Density   1051.94
    INFO> Iteration 3 current conversion 0.405 or 162 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 4 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 4 will generate 38 new bonds
    INFO> Cure_relax initiated on 38 distances (max 0.471 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.422                  0.652
    INFO>          2              0.333                  0.568
    INFO>          3              0.287                  0.534
    INFO>          4              0.222                  0.467
    INFO>          5              0.159                  0.406
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.208 nm x 5.208 nm x 5.208 nm
    INFO> Density                      1052.42
    INFO> Running-average-Density      1032.38
    INFO> Rolling-10-average-Density   1049.88
    INFO> Iteration 4 current conversion 0.500 or 200 bonds

We see here that four CURE iterations were required.  Now comes capping of the 200 oxirane rings that did not react:

.. code-block:: console

    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 200 new bonds
    INFO> Cap_relax initiated on 200 distances (max 0.257 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.199                  0.289
    INFO>          2              0.183                  0.283
    INFO>          3              0.156                  0.286
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.167 nm x 5.167 nm x 5.167 nm
    INFO> Density                      1073.04
    INFO> Running-average-Density      1058.52
    INFO> Rolling-10-average-Density   1074.08
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********

And finally, postcure equilibration and annealing:

.. code-block:: console

    INFO> *************** Postcure in proj-0/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.137 nm x 5.137 nm x 5.137 nm
    INFO> Density                      1091.88
    INFO> Running-average-Density      1087.77
    INFO> Rolling-10-average-Density   1090.81
    INFO> *********** Final data to proj-0/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

So we see our 50% cured system reached about 1.09 g/cc.  I am not claiming any of these equilibration runs are long enough or the system is big enough, but this example just serves to illustrate how ``HTPolyNet`` works.

If we consider the ``hi``-cure system, we see that it took 40 iterations to reach 95% cure.  Below is the end of the console output for that run:

.. code-block:: console

    INFO> Iteration 40 current conversion 0.950 or 380 bonds
    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 20 new bonds
    INFO> Cap_relax initiated on 20 distances (max 0.251 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.191                  0.280
    INFO>          2              0.177                  0.278
    INFO>          3              0.156                  0.271
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.177 nm x 5.177 nm x 5.177 nm
    INFO> Density                      1066.24
    INFO> Running-average-Density      1050.57
    INFO> Rolling-10-average-Density   1065.18
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********
    INFO> *************** Postcure in proj-1/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.147 nm x 5.147 nm x 5.147 nm
    INFO> Density                      1085.58
    INFO> Running-average-Density      1082.20
    INFO> Rolling-10-average-Density   1083.97
    INFO> *********** Final data to proj-1/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************



