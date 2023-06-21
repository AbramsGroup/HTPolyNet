.. _dgeba_results:

Results
-------

The :ref:`results section of the PMS tutorial <pms_results>` does a good job walking you through the project directory structure for any ``HTPolyNet`` run, so we'll skip that here.  

First, if we look in ``proj-0/plots``, we can see the density vs. time of the densification stage:

.. figure:: pics/densification-density.png 

    Density vs. time during the densification stage of the DGEBA/PACM liquid.

As in the polymethylstyrene tutorial, we can make plots of the conversion vs. run time and the cure iteration vs. run time:

.. code-block:: console

    $ htpolynet plots diag --diags diagnostics.log

.. figure:: pics/cure_info.png

    (Left) Conversion vs. wall-clock time; (right) Iteration number vs wall-clock time.

This build was run on a standalone workstation with a very outdated GPU, so it is not surprising it took almost four hours.  The behavior is typical of a build like this:  almost 75% cure is achieved in the first hour of run time, and then nearly *three more hours* are needed to get that last 20%.

Let's generate traces of temperature, density, and potential energy vs time for the entire course of the system build:

.. code-block:: console

    $ htpolynet plots build --proj proj-0 --buildplot t --traces t d p
    
.. figure:: pics/buildtraces.png 

    (Top) Temperature vs. time for 95% cure of DGEBA/PACM; (Middle) Density vs time.  The two topmost plots also report number of bonds; (Bottom) Potential energy vs. time.

Below we show before and after pictures of the system, where all bonds associated with crosslink sites are rendered in "licorice" with the rest in lines, and DGEBA molecules are mauve while PACMs are green:

.. list-table:: 

    * - .. figure:: pics/dge-pac-liq.png

           System before cure.

      - .. figure:: pics/dge-pac-cured.png

           System after cure.