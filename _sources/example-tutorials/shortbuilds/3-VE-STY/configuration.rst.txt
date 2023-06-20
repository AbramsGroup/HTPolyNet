.. _ve_configuration_file:

The Configuration File
======================

We just presented the reactions in the configuration file.  The rest of the file is essentially the same as the original config file for the poly(methyl styrene) system, with the exception of the ``contituents``  and ``densification`` directives:

.. code-block:: yaml

  constituents: {
    STY: { count: 150 },
    GMA: { count: 75 },
    HIE: { stereocenters: [C1, C3] }
  }

Importantly, we can specify GMA here, even though we are not providing an input ``GMA.mol2`` file because GMA is a product that ``HTPolyNet`` recognizes that it will generate by itself.  Secondly, we know that there are two chiral carbons in HIE, which will result in a grand total of 16 distinct diastereromers of GMA being used to build the initial racemic system.

.. code-block:: yaml

  densification: {
    initial_density: 100.0,  # kg/m3
    equilibration: [
      { ensemble: min },
      { ensemble: nvt, temperature: 300, ps: 10.0 },
      { ensemble: npt, temperature: 300, pressure: 10, ps: 100.0, repeat: 8 }
    ]
  }

In the densification, we are stipulating a very low initial density, which is arrived at by trial-and-error because ``gmx insert-molecules`` has problems reliably inserting all the required molecular instances at higher densities, presumably because of GMA's rather long aspect ratio.  To fully densify, we specify that the NPT simulation be repeated 8 times in series, with each step taking 100 ps.  This ensures that gromacs will not complain about the box size changing too much during the densification.  (There may be a better way to do this in a single mdrun execution, but this works.)

Now we are ready to :ref:`run the build <ve_run>`.