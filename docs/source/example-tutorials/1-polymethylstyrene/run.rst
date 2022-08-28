.. _pms_run:

Running the Build
=================

Now, in our working directory ``my_pms/2-polymethylstyrene/``, if you haven't already, go ahead and launch via ``./run.sh``.

.. code-block:: console

    $ ./run.sh

As the two ``htpolynet run`` commands indicate, standard output is being redirected from the console to the file ``lo.log`` for the low-cure run and ``hi.log`` for the high-cure run.  Since the low-cure run executes first, and both use the same ``./lib``, the low-cure run will parameterize and the high-cure run will use those parameterizations.

Console output
^^^^^^^^^^^^^^

In this section, let's go section by section through the console log ``hi.log`` to understand what ``HTPolyNet`` is doing.  The first section just shows a banner and indicates where the library is and where the project directory is, then reports on the available Ambertools software, finally ending by stating the name of the configuration file.

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
    INFO> User library is /home/cfa/htpolynet-tests/test-4/2-polymethylstyrene/lib
    INFO> Starting new build in /home/cfa/htpolynet-tests/test-4/2-polymethylstyrene/proj-1
    INFO> *************************** Ambertools: ***************************
    INFO>   antechamber (ver. 22.0) at antechamber                                        
    INFO>         tleap (ver. 22.0) at tleap                                              
    INFO>      parmchk2 (ver. 22.0) at parmchk2                                           
    INFO> Configuration: pMSTY-hi.yaml

Next comes the monomer and oligomer template parameterizations.  ``HTPolyNet`` first handles the molecules that are explicitly named in the configuration.  Only ``EMB.mol2`` exists; the other two are generated automatically using their respective ``reaction`` directives.  The three molecules implied by the concept of "chaining" are also generated and parameterized.

.. code-block:: console

    INFO> *********** Templates in proj-1/molecules/parameterized ***********
    INFO> 3 molecules explicit in pMSTY-hi.yaml
    INFO> 3 molecules implied by chaining
    INFO> ['EMB~C1=C2~EMB~C1-C2~EMB', 'EMB~C1-C2~EMB~C1=C2~EMB', 'EMB~C1-C2~EMB~C1=C2~EMB~C1-C2~EMB']
    INFO> Generated 6 molecule templates
    INFO> Initial composition is EMB 100
    INFO> 100% conversion is 100 bonds

If we look in ``proj-1/molecules/parameterized`` we'll see the ``gro``, ``itp``, ``top`` and ``grx`` files for each molecule.  (The full parameterizations here were done in the low-cure run since that ``htpolynet`` invocation was first in ``run.sh``.) The first three are Gromacs-specific.  The ``grx`` file contains "extended attributes" of each atom that ``HTPolyNet`` uses internally and are **not** needed for Gromacs.

Next comes initialization of the system topology and coordinates.  Here, using the ``constituents`` directive, ``HTPolyNet`` generates a full system topology and simulation box. The box is filled according to the ``initial_density`` subdirective of the ``densification`` directive.

.. code-block:: console

    INFO> ************** Initialization in proj-1/systems/init **************
    INFO> Topology "init.top" in proj-1/systems/init
    INFO> Initial density: 300.0 kg/m^3
    INFO> Total mass: 1.996e-23 kg
    INFO> Box aspect ratio: 1.0 x 1.0 x 1.0
    INFO> Initial box side lengths: 4.052 nm x 4.052 nm x 4.052 nm
    INFO> Coordinates "init.gro" in proj-1/systems/init
    INFO> Extended attributes "init.grx" in proj-1/systems/init

Next comes a report of the densification of the system.

.. code-block:: console

    INFO> ********** Densification in proj-1/systems/densification **********
    INFO> Running Gromacs: minimization
    INFO> Running Gromacs: nvt ensemble;   10 ps, 300 K
    INFO> Running Gromacs: npt ensemble;  200 ps, 300 K, 10 bar
    INFO> Current box side lengths: 2.872 nm x 2.872 nm x 2.872 nm
    INFO> Density                      842.83
    INFO> Running-average-Density      722.53
    INFO> Rolling-10-average-Density   829.52
    INFO> Densified coordinates in proj-1/systems/densification/densified-npt.gro

Next comes a report of the precure:

.. code-block:: console

    INFO> **************** Precure in proj-1/systems/precure ****************
    INFO> Running Gromacs: npt ensemble;  200 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.892 nm x 2.892 nm x 2.892 nm
    INFO> Density                      825.51
    INFO> Running-average-Density      830.64
    INFO> Rolling-10-average-Density   832.66
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.887 nm x 2.887 nm x 2.887 nm
    INFO> Density                      829.76
    INFO> Running-average-Density      833.70
    INFO> Rolling-10-average-Density   830.74

Next we begin the CURE iterations:

.. code-block:: console

    INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> Attempting to form 95 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 1 will generate 28 new bonds
    INFO> Cure_relax initiated on 28 distances (max 0.455 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.394                  0.615
    INFO>          2              0.341                  0.553
    INFO>          3              0.278                  0.510
    INFO>          4              0.218                  0.451
    INFO>          5              0.162                  0.410
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.820 nm x 2.820 nm x 2.820 nm
    INFO> Density                      886.17
    INFO> Running-average-Density      845.49
    INFO> Rolling-10-average-Density   886.27
    INFO> Iteration 1 current conversion 0.280 or 28 bonds

This first iteration shows that, with a search radius of 0.5 nm, ``HTPolyNet`` identified 28 allowable bonds.  It then forms them and progresses through the relaxation stages until they are at their correct lengths.  Finally it runs the post-iteration NPT MD equilibration, reporting the resulting box dimensions and density.

Next we proceed through CURE iterations 2 through 9:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 19 new bonds
    INFO> Cure_relax initiated on 19 distances (max 0.459 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.400                  0.618
    INFO>          2              0.339                  0.562
    INFO>          3              0.276                  0.511
    INFO>          4              0.219                  0.462
    INFO>          5              0.163                  0.416
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.770 nm x 2.770 nm x 2.770 nm
    INFO> Density                      931.21
    INFO> Running-average-Density      880.93
    INFO> Rolling-10-average-Density   929.96
    INFO> Iteration 2 current conversion 0.470 or 47 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 3 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 3 will generate 20 new bonds
    INFO> Cure_relax initiated on 20 distances (max 0.487 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.436                  0.660
    INFO>          2              0.373                  0.591
    INFO>          3              0.325                  0.549
    INFO>          4              0.262                  0.503
    INFO>          5              0.219                  0.442
    INFO>          6              0.169                  0.421
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.759 nm x 2.759 nm x 2.759 nm
    INFO> Density                      939.34
    INFO> Running-average-Density      896.79
    INFO> Rolling-10-average-Density   938.68
    INFO> Iteration 3 current conversion 0.670 or 67 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 4 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 4 will generate 12 new bonds
    INFO> Cure_relax initiated on 12 distances (max 0.468 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.417                  0.637
    INFO>          2              0.347                  0.564
    INFO>          3              0.285                  0.526
    INFO>          4              0.224                  0.456
    INFO>          5              0.160                  0.417
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.757 nm x 2.757 nm x 2.757 nm
    INFO> Density                      939.41
    INFO> Running-average-Density      900.25
    INFO> Rolling-10-average-Density   941.86
    INFO> Iteration 4 current conversion 0.790 or 79 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 5 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 5 will generate 6 new bonds
    INFO> Cure_relax initiated on 6 distances (max 0.493 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.438                  0.650
    INFO>          2              0.374                  0.604
    INFO>          3              0.323                  0.550
    INFO>          4              0.272                  0.507
    INFO>          5              0.218                  0.457
    INFO>          6              0.159                  0.420
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.766 nm x 2.766 nm x 2.766 nm
    INFO> Density                      929.36
    INFO> Running-average-Density      913.26
    INFO> Rolling-10-average-Density   926.27
    INFO> Iteration 5 current conversion 0.850 or 85 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 6 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 6 will generate 2 new bonds
    INFO> Cure_relax initiated on 2 distances (max 0.444 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.377                  0.614
    INFO>          2              0.321                  0.547
    INFO>          3              0.268                  0.509
    INFO>          4              0.209                  0.436
    INFO>          5              0.158                  0.392
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.737 nm x 2.737 nm x 2.737 nm
    INFO> Density                      958.93
    INFO> Running-average-Density      925.68
    INFO> Rolling-10-average-Density   953.58
    INFO> Iteration 6 current conversion 0.870 or 87 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 7 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 7 will generate 2 new bonds
    INFO> Cure_relax initiated on 2 distances (max 0.479 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.412                  0.620
    INFO>          2              0.355                  0.563
    INFO>          3              0.279                  0.489
    INFO>          4              0.220                  0.454
    INFO>          5              0.153                  0.417
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.758 nm x 2.758 nm x 2.758 nm
    INFO> Density                      937.40
    INFO> Running-average-Density      925.78
    INFO> Rolling-10-average-Density   944.42
    INFO> Iteration 7 current conversion 0.890 or 89 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 8 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Iteration 8 will generate 4 new bonds
    INFO> Cure_drag initiated on 4 distances (max 0.698 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.654
    INFO>          2              0.585
    INFO>          3              0.548
    INFO>          4              0.495
    INFO>          5              0.449
    INFO>          6              0.404
    INFO>          7              0.355
    INFO>          8              0.306
    INFO> Cure_relax initiated on 4 distances (max 0.306 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.265                  0.508
    INFO>          2              0.203                  0.466
    INFO>          3              0.169                  0.420
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.765 nm x 2.765 nm x 2.765 nm
    INFO> Density                      929.01
    INFO> Running-average-Density      886.46
    INFO> Rolling-10-average-Density   928.51
    INFO> Iteration 8 current conversion 0.930 or 93 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 9 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Radius increased to 0.75 nm
    INFO> Iteration 9 will generate 2 new bonds
    INFO> Cure_drag initiated on 2 distances (max 0.633 nm)
    INFO>      Stage  Max-distance (nm)
    INFO>          1              0.587
    INFO>          2              0.549
    INFO>          3              0.500
    INFO>          4              0.438
    INFO>          5              0.403
    INFO>          6              0.345
    INFO>          7              0.307
    INFO> Cure_relax initiated on 2 distances (max 0.307 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.253                  0.484
    INFO>          2              0.211                  0.458
    INFO>          3              0.154                  0.400
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.756 nm x 2.756 nm x 2.756 nm
    INFO> Density                      938.69
    INFO> Running-average-Density      908.04
    INFO> Rolling-10-average-Density   946.97
    INFO> Iteration 9 current conversion 0.950 or 95 bonds

This meets our desired cure of 95\%, so now ``HTPolyNet`` proceeds to capping, 
and not finding any cappable bonds, proceeds to the postcure:

.. code-block:: console

    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 0 new bonds
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********
    INFO> *************** Postcure in proj-1/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.745 nm x 2.745 nm x 2.745 nm
    INFO> Density                      949.71
    INFO> Running-average-Density      946.19
    INFO> Rolling-10-average-Density   949.69
    INFO> *********** Final data to proj-1/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

This just tells us the final density and where the final results are found.  If we look there, we see:

.. code-block:: console

    $ ls -l proj-1/systems/final-results
    final.gro  final.grx  final.top

Now, with the ``gro`` and ``top`` file, you can run whatever Gromacs simulation you like with this system.

