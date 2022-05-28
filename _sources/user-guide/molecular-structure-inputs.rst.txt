Molecular Structure Inputs
~~~~~~~~~~~~~~~~~~~~~~~~~~

HTPolyNet uses the distinct terms "monomer" and "molecule".  A monomer is akin to a "residue" for biomolecules; it has a unique name, and a "molecule" is made of a sequence of one or more monomers.  (Here, "sequence" only refers to the order in which monomers appear in the list of all atoms; not their topological sequence.)  HTPolyNet will build a system out of molecules you specify as input; typically, these would be monomers (i.e., molecules with a single monomer).

Regardless, any molecule you wish to use as a **monomer** must have a mol2 file ``<NAME>.mol2``, where ``<NAME>`` is replaced with the name of the monomer. This name is important; it is how the monomer is called forever inside HTPolyNet.  Input ``*.mol2`` files are all processed via antechamber as part of the parameterization algorithm.

Sample mol2 files for a few monomers are provide in the Library subpackage: ``Library/molecules/inputs``.  You will likely want to create your own.  Most chemical structure drawing programs will output mol2 files.  OpenBabel can also generate them from SMILES strings.  As an example, consider 4,4-diaminodicyclohexylmethane, referred to colloquially as PACM ("pack-em").  The SMILES string for PACM is `C1CC(CCC1CC2CCC(CC2)N)N <https://en.wikipedia.org/wiki/4,4-Diaminodicyclohexylmethane>`_, and using obabel, we can generate a structure for the "PAC" monomer (three-letter monomer names are customary, though not required; they need only be unique identifiers):

.. code-block:: console

    $ echo "C1CC(CCC1CC2CCC(CC2)N)N" | obabel -ismi -h --gen3d -omol2 --title "PAC" | sed s/"UNL1   "/"PAC    "/ > PAC.mol2

The mol2 format HTPolyNet requires for a monomer is minimal and requires only three sections: ``@<TRIPOS>MOLECULE``, ``@<TRIPOS>ATOM`` and ``@<TRIPOS>BOND``.  No matter how you generate a mol2 file, if it corresponds to a molecule that can react, you must edit the mol2 file to give **unique atom names** to any reactive atoms or chiral carbons.  You will then use those unique names to refer specifically to those atoms.

The best way to identify these atoms is by inspecting the 3d structure of the molecule.  Here is the 3d structure of PACM the above command generated:

.. image:: PAC-SC-0-0.png
    :scale: 50 %
    :alt: PACM

Since the two nitrogens are reactive, we will name one "N1" and the other "N2".  It does not matter what they are called so long as the names are unique.  The carbon atoms at the 4 positions of each cyclohexane ring are chiral, so we will name the one connected to N1 "C1" and the other "C2".  These names are changed by directly editing the mol2 file::

    @<TRIPOS>MOLECULE
    PAC
    41 42 0 0 0
    SMALL
    GASTEIGER

    @<TRIPOS>ATOM
         1 C           1.0203    1.1686   -0.4045 C.3     1  PAC        -0.0488
         2 C          -0.3868    1.4530    0.1332 C.3     1  PAC        -0.0375
         3 C2         -0.4239    1.5867    1.6509 C.3     1  PAC         0.0049
         4 C           0.2189    0.3673    2.3129 C.3     1  PAC        -0.0375
         5 C           1.6627    0.1840    1.8377 C.3     1  PAC        -0.0488
         6 C           1.7559    0.0170    0.3181 C.3     1  PAC        -0.0407
         7 C           3.2445   -0.0611   -0.1651 C.3     1  PAC        -0.0474
         8 C           4.0849   -1.2509    0.3999 C.3     1  PAC        -0.0407
         9 C           5.5341   -1.2664   -0.1535 C.3     1  PAC        -0.0488
        10 C           6.3098   -2.5522    0.1636 C.3     1  PAC        -0.0375
        11 C1          5.4974   -3.8029   -0.1700 C.3     1  PAC         0.0049
        12 C           4.1937   -3.8000    0.6212 C.3     1  PAC        -0.0375
        13 C           3.3524   -2.5924    0.2247 C.3     1  PAC        -0.0488
        14 N1          6.2599   -5.0172    0.1162 N.3     1  PAC        -0.3272
        15 N2         -1.8168    1.7369    2.0786 N.3     1  PAC        -0.3272
        16 H           1.6047    2.0898   -0.3424 H       1  PAC         0.0268
        17 H           0.9019    0.9202   -1.4627 H       1  PAC         0.0268
        18 H          -1.0564    0.6483   -0.1927 H       1  PAC         0.0280
        19 H          -0.7633    2.3773   -0.3343 H       1  PAC         0.0280
        20 H           0.1247    2.4885    1.9532 H       1  PAC         0.0458
        21 H          -0.3534   -0.5388    2.0744 H       1  PAC         0.0280
        22 H           0.2067    0.4761    3.4022 H       1  PAC         0.0280
        23 H           2.0776   -0.7078    2.3325 H       1  PAC         0.0268
        24 H           2.2691    1.0366    2.1678 H       1  PAC         0.0268
        25 H           1.2371   -0.9012    0.0434 H       1  PAC         0.0301
        26 H           3.7432    0.8605    0.1294 H       1  PAC         0.0271
        27 H           3.2593   -0.0975   -1.2596 H       1  PAC         0.0271
        28 H           4.1835   -1.0879    1.4730 H       1  PAC         0.0301
        29 H           6.0813   -0.4176    0.2686 H       1  PAC         0.0268
        30 H           5.5482   -1.1352   -1.2427 H       1  PAC         0.0268
        31 H           6.5982   -2.5580    1.2292 H       1  PAC         0.0280
        32 H           7.2463   -2.5515   -0.4099 H       1  PAC         0.0280
        33 H           5.2588   -3.8065   -1.2451 H       1  PAC         0.0458
        34 H           4.3975   -3.7699    1.7048 H       1  PAC         0.0280
        35 H           3.6232   -4.7131    0.4347 H       1  PAC         0.0280
        36 H           2.4417   -2.5989    0.8401 H       1  PAC         0.0268
        37 H           3.0264   -2.7222   -0.8175 H       1  PAC         0.0268
        38 H           6.5386   -5.0249    1.0976 H       1  PAC         0.1185
        39 H           7.1205   -5.0120   -0.4231 H       1  PAC         0.1185
        40 H          -2.3522    0.9246    1.7729 H       1  PAC         0.1185
        41 H          -2.2309    2.5311    1.5900 H       1  PAC         0.1185
    @<TRIPOS>BOND
         1     1     2    1
         2     2     3    1
         3     3     4    1
         4     4     5    1
         5     5     6    1
         6     1     6    1
         7     6     7    1
         8     7     8    1
         9     8     9    1
        10     9    10    1
        11    10    11    1
        12    11    12    1
        13    12    13    1
        14     8    13    1
        15    11    14    1
        16     3    15    1
        17     1    16    1
        18     1    17    1
        19     2    18    1
        20     2    19    1
        21     3    20    1
        22     4    21    1
        23     4    22    1
        24     5    23    1
        25     5    24    1
        26     6    25    1
        27     7    26    1
        28     7    27    1
        29     8    28    1
        30     9    29    1
        31     9    30    1
        32    10    31    1
        33    10    32    1
        34    11    33    1
        35    12    34    1
        36    12    35    1
        37    13    36    1
        38    13    37    1
        39    14    38    1
        40    14    39    1
        41    15    40    1
        42    15    41    1
