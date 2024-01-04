.. _dgeba_run:

Running the Build
=================

Now, in our working directory ``my_dgeba_pacm/4-pacm-dgeba-epoxy-thermoset/``, if you haven't already, go ahead and launch via ``./run.sh``.

.. code-block:: console

    $ ./run.sh

As the two ``htpolynet run`` commands indicate, standard output is being redirected from the console to the file ``console.log``.

First the banner and the runtime initialization messages:

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
    INFO> User library is /home/cfa/htpolynet-tutorials/1.0.8/4-pacm-dgeba-epoxy-thermoset/lib
    INFO> New project in /home/cfa/htpolynet-tutorials/1.0.8/4-pacm-dgeba-epoxy-thermoset/proj-0
    INFO> *************************** Ambertools: ***************************
    INFO> ********************  antechamber (ver. 22.0) *********************
    INFO> ********************        tleap (ver. 22.0) *********************
    INFO> ********************     parmchk2 (ver. 22.0) *********************
    INFO> Configuration: DGEPAC.yaml

Then a description of all the molecules that need to be created:

.. code-block:: console

    INFO> *********** Templates in proj-0/molecules/parameterized ***********
    INFO> 22 molecules detected in DGEPAC.yaml
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

Notice that the symmetry-equivalent atoms in both monomers result in a total of 12 unique cure reactions forming secondary and tertiary amines, as well as two capping reactions for DGEBA.  Now the initial construction of the system and its densification are next:

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
    INFO> Current box side lengths: 5.225 nm x 5.225 nm x 5.225 nm
    INFO> Density                      1046.77
    INFO> Running-average-Density       875.31
    INFO> Rolling-10-average-Density   1044.65
    INFO> Densified coordinates in proj-0/systems/densification/densified-npt.gro
    INFO> **************** Precure in proj-0/systems/precure ****************
    INFO> Running Gromacs: npt ensemble; 200.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.219 nm x 5.219 nm x 5.219 nm
    INFO> Density                      1050.62
    INFO> Running-average-Density      1047.96
    INFO> Rolling-10-average-Density   1050.99
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.224 nm x 5.224 nm x 5.224 nm
    INFO> Density                      1047.39
    INFO> Running-average-Density      1048.71
    INFO> Rolling-10-average-Density   1049.57

Note that we achieve a pretty good initial density for this liquid of about 1.05 g/cc.  Now the CURE algorithm begins, seeking to form 380 out of the total possible 400 bonds:

.. code-block:: console

    INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> Attempting to form 380 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    
The run then enters the CURE iterations.  The console output shows that 121 bonds form in just those first two iterations:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 1 will generate 71 new bonds
    INFO> Step "cure_relax" initiated on 71 distances (max 0.478 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.418                  0.639
    INFO>          2              0.349                  0.588
    INFO>          3              0.279                  0.526
    INFO>          4              0.224                  0.461
    INFO>          5              0.156                  0.405
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.228 nm x 5.228 nm x 5.228 nm
    INFO> Density                      1043.18
    INFO> Running-average-Density      1022.51
    INFO> Rolling-10-average-Density   1042.26
    INFO> Iteration 1 current conversion 0.177 or 71 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 50 new bonds
    INFO> Step "cure_relax" initiated on 50 distances (max 0.489 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.448                  0.662
    INFO>          2              0.362                  0.586
    INFO>          3              0.325                  0.546
    INFO>          4              0.257                  0.487
    INFO>          5              0.212                  0.446
    INFO>          6              0.157                  0.403
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.216 nm x 5.216 nm x 5.216 nm
    INFO> Density                      1049.56
    INFO> Running-average-Density      1025.18
    INFO> Rolling-10-average-Density   1045.13
    INFO> Iteration 2 current conversion 0.302 or 121 bonds

This rate of bond formation doesn't last, though.  It gets harder and harder to find potential partners the deeper into the cure we go.  By the end, this build required 32 CURE iterations.  The console output for the last two is below, and we see both required pre-bond dragging before bond formation and relaxation:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 31 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 31 will generate 5 new bonds
    INFO> Step "cure_drag" initiated on 5 distances (max 0.994 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.926
    INFO>          2              0.886
    INFO>          3              0.817
    INFO>          4              0.764
    INFO>          5              0.707
    INFO>          6              0.653
    INFO>          7              0.593
    INFO>          8              0.536
    INFO>          9              0.478
    INFO>         10              0.423
    INFO>         11              0.361
    INFO>         12              0.303
    INFO> Step "cure_relax" initiated on 5 distances (max 0.303 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.261                  0.499
    INFO>          2              0.207                  0.444
    INFO>          3              0.154                  0.400
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.199 nm x 5.199 nm x 5.199 nm
    INFO> Density                      1053.87
    INFO> Running-average-Density      1040.71
    INFO> Rolling-10-average-Density   1055.13
    INFO> Iteration 31 current conversion 0.945 or 378 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 32 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 32 will generate 2 new bonds
    INFO> Step "cure_drag" initiated on 2 distances (max 0.720 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.659
    INFO>          2              0.627
    INFO>          3              0.574
    INFO>          4              0.528
    INFO>          5              0.484
    INFO>          6              0.447
    INFO>          7              0.389
    INFO>          8              0.346
    INFO>          9              0.305
    INFO> Step "cure_relax" initiated on 2 distances (max 0.305 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.256                  0.490
    INFO>          2              0.208                  0.447
    INFO>          3              0.157                  0.406
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.192 nm x 5.192 nm x 5.192 nm
    INFO> Density                      1057.63
    INFO> Running-average-Density      1043.42
    INFO> Rolling-10-average-Density   1056.08
    INFO> Iteration 32 current conversion 0.950 or 380 bonds

With this conversion reached, now comes capping of the 20 oxirane rings that did not react:

.. code-block:: console

    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 20 new bonds
    INFO> Step "cap_relax" initiated on 20 distances (max 0.252 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.189                  0.277
    INFO>          2              0.174                  0.281
    INFO>          3              0.155                  0.279
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.192 nm x 5.192 nm x 5.192 nm
    INFO> Density                      1057.48
    INFO> Running-average-Density      1048.40
    INFO> Rolling-10-average-Density   1058.75
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********

And finally, postcure equilibration and annealing:

.. code-block:: console

    INFO> *************** Postcure in proj-0/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.146 nm x 5.146 nm x 5.146 nm
    INFO> Density                      1085.93
    INFO> Running-average-Density      1085.08
    INFO> Rolling-10-average-Density   1087.24
    INFO> *********** Final data to proj-0/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

So we see our 95% cured system reached a density about 1.09 g/cc.  I am not claiming any of these equilibration runs are long enough or the system is big enough, but this example just serves to illustrate how ``HTPolyNet`` works.

Let's now take a look at some selected :ref:`results <dgeba_results>`.