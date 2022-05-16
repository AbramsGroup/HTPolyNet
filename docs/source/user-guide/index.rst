User Guide
==========

Overview
~~~~~~~~

HTPolyNet generally follows this procedure:

1. If necessary, GAFF-parameterize any input monomers, and build any intermediate molecules dictated by the polymerization chemistry and GAFF-parameterize those.  Outputs of parameterizations are typically saved in a Library so that repeated parameterizations are not done.
2. Based on a specified composition, generate an initial simulation box, and equilibrate it to a liquid-like density.
3. Perform CURE (Connect-Update-Relax-Equilibrate) iterations to introduce intermonomer bonds according to the input chemistry until a desired conversion is met or an execution threshhold is reached.
4. Finalize by performing any requested post-cure chemistry.

An HTPolyNet run therefore requires two types of inputs: molecular structures of monomers, and a configuration file.  These are described below.

Inputs
~~~~~~

HTPolyNet requires minimally two types of Inputs

1. Molecular structures of any monomers as Tripos mol2 files; these are used directly as inputs to ``antechamber``.
2. An input configuration file in ``YAML`` format for controlling the ``htpolynet`` run.

Molecular structures
^^^^^^^^^^^^^^^^^^^^

HTPolyNet uses the distinct terms "monomer" and "molecule".  A monomer is akin to a "residue" for biomolecules; it has a unique name, and a "molecule" is made of a sequence of one or more monomers.  (Here, "sequence" only refers to the order in which monomers appear in the list of all atoms; not their topological sequence.)  HTPolyNet will build a system out of molecules you specify as input; typically, these would be monomers (i.e., molecules with a single monomer).

Regardless, any molecule you wish to use as a monomer must have a mol2 file ``<NAME>.mol2``, where ``<NAME>`` is replaced with the name of the monomer. This name is important; it is how the monomer is called forever inside HTPolyNet.

Sample mol2 files for a few monomers are provide in the Library subpackage: ``Library/molecules/inputs``.  You will likely want to create your own.  Most chemical structure drawing programs will output mol2 files.  OpenBabel can also generate them from SMILES strings; e.g., for butane:

.. code-block:: console

    $  echo "CCCC" | obabel -ismi -omol2 --gen3d --title BUTANE | sed s/"UNL1  "/"BUTANE"/ > BUTANE.mol2

So to generate mol2 files for your monomers, you have a few options.

The mol2 format HTPolyNet requires for a monomer is minimal and requires only three sections: ``@<TRIPOS>MOLECULE``, ``@<TRIPOS>ATOM`` and ``@<TRIPOS>BOND``.  For a molecule with more than one monomer, a ``@<TRIPOS>SUBSTRUCTURE`` section is also required.  Below is the mol2 file for butane we just created above::

    @<TRIPOS>MOLECULE
    BUTANE
    14 13 0 0 0
    SMALL
    GASTEIGER

    @<TRIPOS>ATOM
         1 C           0.9926    0.1046    0.0014 C.3     1  BUTANE     -0.0653
         2 C           2.5129    0.0970   -0.0099 C.3     1  BUTANE     -0.0562
         3 C           3.0598   -0.7349   -1.1682 C.3     1  BUTANE     -0.0562
         4 C           4.5801   -0.7416   -1.1802 C.3     1  BUTANE     -0.0653
         5 H           0.6228    0.7054    0.8382 H       1  BUTANE      0.0230
         6 H           0.5960   -0.9100    0.1102 H       1  BUTANE      0.0230
         7 H           0.5955    0.5315   -0.9250 H       1  BUTANE      0.0230
         8 H           2.8777   -0.3052    0.9424 H       1  BUTANE      0.0263
         9 H           2.8771    1.1281   -0.0871 H       1  BUTANE      0.0263
        10 H           2.6944   -0.3331   -2.1205 H       1  BUTANE      0.0263
        11 H           2.6961   -1.7661   -1.0906 H       1  BUTANE      0.0263
        12 H           4.9760    0.2731   -1.2896 H       1  BUTANE      0.0230
        13 H           4.9499   -1.3426   -2.0168 H       1  BUTANE      0.0230
        14 H           4.9778   -1.1679   -0.2537 H       1  BUTANE      0.0230
    @<TRIPOS>BOND
         1     1     2    1
         2     2     3    1
         3     3     4    1
         4     1     5    1
         5     1     6    1
         6     1     7    1
         7     2     8    1
         8     2     9    1
         9     3    10    1
        10     3    11    1
        11     4    12    1
        12     4    13    1
        13     4    14    1


Configuration files
^^^^^^^^^^^^^^^^^^^

That
