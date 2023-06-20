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
    INFO> User library is /home/cfa/htpolynet-tutorials/v1.0.6/4-pacm-dgeba-epoxy-thermoset/lib
    INFO> New project in /home/cfa/htpolynet-tutorials/v1.0.6/4-pacm-dgeba-epoxy-thermoset/proj-0
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
    INFO> Current box side lengths: 5.245 nm x 5.245 nm x 5.245 nm
    INFO> Density                      1034.74
    INFO> Running-average-Density       862.50
    INFO> Rolling-10-average-Density   1037.43
    INFO> Densified coordinates in proj-0/systems/densification/densified-npt.gro
    INFO> **************** Precure in proj-0/systems/precure ****************
    INFO> Running Gromacs: npt ensemble; 200.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.228 nm x 5.228 nm x 5.228 nm
    INFO> Density                      1045.20
    INFO> Running-average-Density      1043.63
    INFO> Rolling-10-average-Density   1045.36
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.219 nm x 5.219 nm x 5.219 nm
    INFO> Density                      1050.51
    INFO> Running-average-Density      1046.98
    INFO> Rolling-10-average-Density   1048.83

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
    INFO> Iteration 1 will generate 74 new bonds
    INFO> Step "cure_relax" initiated on 74 distances (max 0.471 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.399                  0.629
    INFO>          2              0.344                  0.574
    INFO>          3              0.285                  0.508
    INFO>          4              0.219                  0.455
    INFO>          5              0.160                  0.413
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.222 nm x 5.222 nm x 5.222 nm
    INFO> Density                      1047.15
    INFO> Running-average-Density      1022.39
    INFO> Rolling-10-average-Density   1049.11
    INFO> Iteration 1 current conversion 0.185 or 74 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 47 new bonds
    INFO> Step "cure_relax" initiated on 47 distances (max 0.490 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.437                  0.673
    INFO>          2              0.384                  0.612
    INFO>          3              0.318                  0.535
    INFO>          4              0.263                  0.491
    INFO>          5              0.215                  0.453
    INFO>          6              0.156                  0.404
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.209 nm x 5.209 nm x 5.209 nm
    INFO> Density                      1053.69
    INFO> Running-average-Density      1025.92
    INFO> Rolling-10-average-Density   1051.94
    INFO> Iteration 2 current conversion 0.302 or 121 bonds

This rate of bond formation doesn't last, though.  It gets harder and harder to find potential partners the deeper into the cure we go.  By the end, this build required 45 CURE iterations.  The console output for the last two is below, and we see they each generated only one bond, and both required pre-bond dragging before bond formation and relaxation:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 44 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Radius increased to 1.0 nm
    INFO> Iteration 44 will generate 1 new bond
    INFO> Step "cure_drag" initiated on 1 distance (max 0.927 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.887
    INFO>          2              0.813
    INFO>          3              0.752
    INFO>          4              0.700
    INFO>          5              0.647
    INFO>          6              0.594
    INFO>          7              0.524
    INFO>          8              0.479
    INFO>          9              0.415
    INFO>         10              0.360
    INFO>         11              0.299
    INFO> Step "cure_relax" initiated on 1 distance (max 0.299 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.257                  0.480
    INFO>          2              0.197                  0.435
    INFO>          3              0.150                  0.396
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.188 nm x 5.188 nm x 5.188 nm
    INFO> Density                      1060.00
    INFO> Running-average-Density      1042.28
    INFO> Rolling-10-average-Density   1059.91
    INFO> Iteration 44 current conversion 0.948 or 379 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 45 begins ~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Radius increased to 1.0 nm
    INFO> Iteration 45 will generate 1 new bond
    INFO> Step "cure_drag" initiated on 1 distance (max 0.892 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.845
    INFO>          2              0.788
    INFO>          3              0.725
    INFO>          4              0.682
    INFO>          5              0.610
    INFO>          6              0.571
    INFO>          7              0.514
    INFO>          8              0.464
    INFO>          9              0.415
    INFO>         10              0.357
    INFO>         11              0.307
    INFO> Step "cure_relax" initiated on 1 distance (max 0.307 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.264                  0.489
    INFO>          2              0.212                  0.440
    INFO>          3              0.154                  0.385
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.187 nm x 5.187 nm x 5.187 nm
    INFO> Density                      1060.60
    INFO> Running-average-Density      1049.17
    INFO> Rolling-10-average-Density   1061.22
    INFO> Iteration 45 current conversion 0.950 or 380 bonds

With this conversion reached, now comes capping of the 20 oxirane rings that did not react:

.. code-block:: console

    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 20 new bonds
    INFO> Step "cap_relax" initiated on 20 distances (max 0.249 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.194                  0.285
    INFO>          2              0.179                  0.279
    INFO>          3              0.153                  0.281
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.180 nm x 5.180 nm x 5.180 nm
    INFO> Density                      1064.91
    INFO> Running-average-Density      1052.92
    INFO> Rolling-10-average-Density   1064.26
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********

And finally, postcure equilibration and annealing:

.. code-block:: console

    INFO> *************** Postcure in proj-0/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble; 100.00 ps,  300.00 K,   1.00 bar
    INFO> Current box side lengths: 5.151 nm x 5.151 nm x 5.151 nm
    INFO> Density                      1082.81
    INFO> Running-average-Density      1080.81
    INFO> Rolling-10-average-Density   1083.31
    INFO> *********** Final data to proj-0/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

So we see our 95% cured system reached a density about 1.08 g/cc.  I am not claiming any of these equilibration runs are long enough or the system is big enough, but this example just serves to illustrate how ``HTPolyNet`` works.

Let's now take a look at some selected :ref:`results <dgeba_results>`.