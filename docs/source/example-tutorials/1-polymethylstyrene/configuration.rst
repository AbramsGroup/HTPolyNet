.. _pms_configuration_file:

The Configuration File
======================

The example you fetched has two configuration files that are essentially identical, except one sets the desired cure to 95\% (hi) and the other only to 50\% (lo).  We previously described an example `configuration file for making a system of polystyrene <configuration_files>`, where all the directives were explained in detail.  Here we will only highlight *excerpts* from ``pMSTY-hi.yaml``.

Constituents Directives
~~~~~~~~~~~~~~~~~~~~~~~

We'll use a very small system of only 100 monomers.

.. code-block:: yaml

  constituents: {
    EMB: {count: 100}
  }


Reactions Directives
~~~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    reactions:
      - {
          name:        'EMB1_1',
          stage:       cure,
          reactants:   {1: EMB, 2: EMB},
          product:     EMB~C1-C2~EMB,
          probability: 1.0,
          atoms: {
            A: {reactant: 1, resid: 1, atom: C1, z: 1},
            B: {reactant: 2, resid: 1, atom: C2, z: 1}
          },
          bonds: [
            {atoms: [A, B], order: 1}
          ]
        }
      - {
          name:         'EMBCC',
          stage:        post-cure,
          reactants:    {1: EMB},
          product:      EMBCC,
          probability:  1.0,
          atoms: {
            A: {reactant: 1, resid: 1, atom: C1, z: 1},
            B: {reactant: 1, resid: 1, atom: C2, z: 1}
          },
          bonds: [
            {atoms: [A, B], order: 2}
          ]
        }
