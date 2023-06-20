.. _pms_configuration_file:

The Configuration File
======================

We previously described an example :ref:`configuration file for making a system of polystyrene <configuration_files>`, where all the directives were explained in detail.  Here we will only highlight *excerpts* from ``pMSTY.yaml``.

Constituents Directives
~~~~~~~~~~~~~~~~~~~~~~~

We'll use a very small system of only 127 monomers.

.. code-block:: yaml

  constituents:
    EMB:
      count: 127

Reactions Directives
~~~~~~~~~~~~~~~~~~~~

We need only define the polymerization reaction that links two monomers and the capping reaction that "unactivates" any unreacted monomers.

.. code-block:: yaml

  reactions:
    - name:        emb1_1
      stage:       cure
      reactants:
        1: EMB
        2: EMB
      product:     EMB~C1-C2~EMB
      probability: 1.0
      atoms:
          A:
            reactant: 1
            resid: 1
            atom: C1
            z: 1
          B:
            reactant: 2
            resid: 1
            atom: C2
            z: 1
      bonds:
        - atoms:
            - A
            - B
          order: 1
    - name:         embCC
      stage:        cap
      reactants:
        1: EMB
      product:      EMBCC
      probability:  1.0
      atoms:
        A:
          reactant: 1
          resid: 1
          atom: C1
          z: 1
        B:
          reactant: 1
          resid: 1
          atom: C2
          z: 1
      bonds:
        - atoms:
            - A
            - B
          order: 2
