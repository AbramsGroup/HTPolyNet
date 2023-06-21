.. _tutorials_postsim_analyses:

Examples Using ``htpolynet postsim`` and ``analyze``
====================================================

Here, we will return to the DGEBA/PACM build we did in the earlier tutorials.  Navigate into that base directory and create the file ``postsim.yaml`` with these contents:

.. code-block:: yaml

    - anneal:
        input_top: systems/final-results/final.top
        input_gro: systems/final-results/final.gro
        P: 1
        T0: 300
        T1: 600
        ncycles: 2
        T0_to_T1_ps: 10
        T1_ps: 10
        T1_to_T0_ps: 10
        T0_ps: 10
    - equilibrate:
        input_top: systems/final-results/final.top
        input_gro: postsim/anneal/anneal.gro
        T: 300
        ps: 10
    - ladder:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/ladder-heat
        Tlo: 300
        Thi: 600
        deltaT: 5
        ps_per_rise: 10
        ps_per_run: 90
        warmup_ps: 10 
    - ladder:
        input_top: systems/final-results/final.top
        input_gro: postsim/ladder-heat/ladder.gro
        subdir: postsim/ladder-cool
        Tlo: 300
        Thi: 600
        deltaT: -5
        ps_per_rise: 10
        ps_per_run: 90
        warmup_ps: 10 
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-x
        T: 300
        P: 1
        direction: x
        edot: 0.001
        ps: 10
        gromacs:
          gmx: gmx
          mdrun: gmx mdrun -ntmpi 1
          options: -quiet -nobackup
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-y
        T: 300
        P: 1
        direction: y
        edot: 0.001
        ps: 10
        gromacs:
          gmx: gmx
          mdrun: gmx mdrun -ntmpi 1
          options: -quiet -nobackup
    - deform:
        input_top: systems/final-results/final.top
        input_gro: postsim/equilibrate/equilibrate.gro
        subdir: postsim/deform-z
        T: 300
        P: 1
        direction: z
        edot: 0.001
        ps: 10
        gromacs:
          gmx: gmx
          mdrun: gmx mdrun -ntmpi 1
          options: -quiet -nobackup

This file will instruct HTPolyNet to use Gromacs to conduct several MD simulations.  First, the final results of the build are annealed for two cycles between 300 and 600 K.  Then the result of the annealing is equilibrated at 300 K.  Following, both an increasing and then decreasing temperature ladders are run; in parallel (using the equilibrated output), three uniaxial deformations are run.

*Important disclaimer*:  Although ``htpolynet postsim`` is meant to facilitate production-grade MD simulations, those illustrated here are NOT; they are WAY TOO SHORT.  It is your responsibility to determine how long a simulation is needed for each kind of simulation.

(in progress)
