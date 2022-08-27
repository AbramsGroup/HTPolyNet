.. _dgeba_configuration_file:

The Configuration File
======================

The example you fetched has two configuration files that are essentially identical, except one sets the desired cure to 95\% (hi) and the other only to 50\% (lo).  We previously described an example :ref:`configuration file for making a system of polystyrene <configuration_files>`, where all the directives were explained in detail.  Here we will only highlight the ``constituents`` directive.

We'll use a very small, stoichiometric system of only 200 DGEBA's and 100 PACM's.

.. code-block:: yaml

  constituents: {  
    DGE: { count: 200, symmetry_equivalent_atoms: [[C1,C2],[C3,C4],[O1,O2]], stereocenters: [C3] },                                       
    PAC: { count: 100, symmetry_equivalent_atoms: [[N1,N2],[C1,C2]], stereocenters: [C1] }                                               
  }

Notice that the ``contituents`` directive is where we specify both symmetry-equivalent atoms and stereocenters in each monomer by referencing atom names.  DGEBA has three sets of symmetry-equivalent atoms and one symmetry-unique stereocenter, for example.
