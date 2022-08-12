.. _pms_run:

Running the Build
=================

Now, in our working directory ``my_pms``, if you haven't already, go ahead and launch via ``./run.sh``.

.. code-block:: console

    $ ./run.sh

As the two ``htpolynet run`` commands indicate, standard output is being redirected from the console to the file ``lo.log`` for the low-cure run and ``hi.log`` for the high-cure run.  Let's consider the contents of ``lo.log`` and then examine the project directory that was created.

Console output
^^^^^^^^^^^^^^

In this section, let's go section by section through the console log to understand what ``HTPolyNet`` is doing.  The first section just shows a banner and indicates where the library is and where the project directory is, then reports on the available Ambertools software, finally ending by stating the name of the configuration file.

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
    INFO> User library is /home/cfa/htpolynet-tests/test-2/2-polymethylstyrene/lib
    INFO> Starting new build in /home/cfa/htpolynet-tests/test-2/2-polymethylstyrene/proj-0
    INFO> *************************** Ambertools: ***************************
    INFO>   antechamber (ver. 22.0) at antechamber                                        
    INFO>         tleap (ver. 22.0) at tleap                                              
    INFO>      parmchk2 (ver. 22.0) at parmchk2                                           
    INFO> Configuration: pMSTY-lo.yaml

Next comes the monomer and oligomer template parameterizations.  ``HTPolyNet`` first handles the molecules that are explicitly named in the configuration.  Only ``EMB.mol2`` exists; the other two are generated automatically using their respective ``reaction`` directives.  The three molecules implied by the concept of "chaining" are also generated and parameterized.

.. code-block:: console

    INFO> *********** Templates in proj-0/molecules/parameterized ***********
    INFO> 3 molecules explicit in pMSTY-lo.yaml
    INFO> AmberTools> parameterizing EMB.mol2
    INFO> AmberTools> parameterizing EMB~C1-C2~EMB.mol2
    INFO> AmberTools> parameterizing EMBCC.mol2
    INFO> 3 molecules implied by chaining
    INFO> ['EMB~C1=C2~EMB~C1-C2~EMB', 'EMB~C1-C2~EMB~C1=C2~EMB', 'EMB~C1-C2~EMB~C1=C2~EMB~C1-C2~EMB']
    INFO> AmberTools> parameterizing EMB~C1=C2~EMB~C1-C2~EMB.mol2
    INFO> AmberTools> parameterizing EMB~C1-C2~EMB~C1=C2~EMB.mol2
    INFO> AmberTools> parameterizing EMB~C1-C2~EMB~C1=C2~EMB~C1-C2~EMB.mol2
    INFO> Generated 6 molecule templates
    INFO> Initial composition is EMB 100
    INFO> 100% conversion is 100 bonds

If we look in ``proj-0/molecules/parameterized`` we'll see a LOT of files.  What's important are the ``gro``, ``itp``, ``top`` and ``grx`` files for each molecule.  The first three are Gromacs-specific.  The ``grx`` file contains "extended attributes" of each atom that ``HTPolyNet`` uses internally and are NOT needed for Gromacs.

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
    INFO> Running Gromacs: nvt ensemble;   10 ps, 300 K
    INFO> Running Gromacs: npt ensemble;  200 ps, 300 K, 10 bar
    INFO> Current box side lengths: 2.890 nm x 2.890 nm x 2.890 nm
    INFO> Density                       826.996948
    INFO> Running-average-Density       728.728314
    INFO> Rolling-10-average-Density    828.510007
    INFO> Densified coordinates in proj-0/systems/densification/densified-npt.gro

Next comes a report of the precure:

.. code-block:: console

    INFO> **************** Precure in proj-0/systems/precure ****************
    INFO> Running Gromacs: npt ensemble;  200 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.887 nm x 2.887 nm x 2.887 nm
    INFO> Density                       829.116699
    INFO> Running-average-Density       828.608087
    INFO> Rolling-10-average-Density    831.531906
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.885 nm x 2.885 nm x 2.885 nm
    INFO> Density                       830.845398
    INFO> Running-average-Density       833.242340
    INFO> Rolling-10-average-Density    831.299469

Next we begin the CURE iterations:

.. code-block:: console

    INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********
    INFO> Attempting to form 50 bonds
    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 1 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 1 will generate 31 new bonds
    INFO> Cure_relax initiated on 31 distances (max 0.485 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.428                  0.639
    INFO>          2              0.397                  0.620
    INFO>          3              0.314                  0.538
    INFO>          4              0.263                  0.509
    INFO>          5              0.212                  0.444
    INFO>          6              0.163                  0.420
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.814 nm x 2.814 nm x 2.814 nm
    INFO> Density                       891.040100
    INFO> Running-average-Density       837.843547
    INFO> Rolling-10-average-Density    887.186151
    INFO> Iteration 1 current conversion 0.310 or 31 bonds

This first iteration shows that, with a search radius of 0.5 nm, ``HTPolyNet`` identified 31 allowable bonds that could form.  It then forms them and progresses through the relaxation stages until they are at their correct lenghts.  Finally it runs the post-iteration NPT MD equilibration, reporting the resulting box dimensions and density.

Next we proceed through a second CURE iteration:

.. code-block:: console

    INFO> ~~~~~~~~~~~~~~~~~~~~~~~ Iteration 2 begins ~~~~~~~~~~~~~~~~~~~~~~~~
    INFO> Bond search using radius 0.5 nm initiated
    INFO> Iteration 2 will generate 20 new bonds
    INFO> Cure_relax initiated on 20 distances (max 0.444 nm)
    INFO>      Stage  Max-distance (nm)  Max-1-4-distance (nm)
    INFO>          1              0.388                  0.624
    INFO>          2              0.340                  0.574
    INFO>          3              0.271                  0.511
    INFO>          4              0.221                  0.464
    INFO>          5              0.163                  0.412
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.770 nm x 2.770 nm x 2.770 nm
    INFO> Density                       931.056519
    INFO> Running-average-Density       890.387364
    INFO> Rolling-10-average-Density    929.248755
    INFO> Iteration 2 current conversion 0.510 or 51 bonds
    INFO> ************************* Capping begins **************************
    INFO> Capping will generate 0 new bonds
    INFO> ********** Connect-Update-Relax-Equilibrate (CURE) ends ***********

Since the conversion exceeds the specified desired conversion of 50\%, ``HTPolyNet`` proceeds to capping.  Here, it found no unreacted styrenes that needed capping, so no new bonds were generated here.  Next, ``HTPolyNet`` proceeds to postcure:

.. code-block:: console

    INFO> *************** Postcure in proj-0/systems/postcure ***************
    INFO> Annealing: 5 points for 2 cycles over 160 ps
    INFO> Annealed coordinates in annealed.gro
    INFO> Running Gromacs: npt ensemble;  100 ps, 300 K, 1 bar
    INFO> Current box side lengths: 2.763 nm x 2.763 nm x 2.763 nm
    INFO> Density                       937.608765
    INFO> Running-average-Density       937.181557
    INFO> Rolling-10-average-Density    929.233405
    INFO> *********** Final data to proj-0/systems/final-results ************
    INFO> ********************* HTPolyNet runtime ends **********************

This just tells us the final density and where the final results are found.  If we look there, we see:

.. code-block:: console

    $ ls -l proj-0/systems/final-results
    final.gro  final.grx  final.top

Now, with the ``gro`` and ``top`` file, you can run whatever Gromacs simulation you like with this system.

Overall behavior
^^^^^^^^^^^^^^^^

Using ``htpolynet plots`` we can generate a few interesting graphics that help characterize a build.  In this tutorial, we generated a low-cure build under ``proj-0`` and a high-cure build under ``proj-1``.  Diagnostic output for each run is in ``diagnostics-lo.log`` and ``diagnostics-hi.log``, respectively.

First, we can make plots of the conversion vs. run time and the cure iteration vs. run time:

.. code-block:: console

    $ htpolynet plots -logs diagnostics-*.log

This generates ``cure-info.png``: 

.. image:: cure-info.png 

We can see here that the 95\% cure took about 8 and a half minutes of run time (which is not really impressive since this is a **very** small system).  Fully two-thirds of the run time is consumed realizing the final 15\% of the cure.

Second, we can make plots that track the temperature and density throughout the entire build process:

.. code-block:: console

    $ htpolynet plots -proj proj-0

