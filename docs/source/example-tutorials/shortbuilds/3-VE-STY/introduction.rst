.. _ve_introduction:

Introduction
============

The first step is to get to a clean working directory, and then use ``htpolynet fetch-example`` to setup this example.

.. code-block:: console

   $ mkdir my_vinyl_ester
   $ cd my_vinyl_ester
   $ htpolynet fetch-example -n 3
   $ ls
   3-bisgma-styrene-thermoset/
   $ tree 3-bisgma-styrene-thermoset
   3-bisgma-styrene-thermoset
   ├── GMASTY.yaml
   ├── lib
   │   └── molecules
   │       ├── inputs
   │       └── parameterized
   ├── README.md
   └── run.sh

The bash script ``run.sh`` is a *suggested* way to run a build, and it contains the commands to generate the input ``mol2`` files for activated bisGMA and styrene.

.. code-block:: bash

   #!/bin/bash -l

   mollib='./lib/molecules'
   conda activate mol-env
   rm -rf proj-* *log ${mollib}/parameterized/*

   STYRENE="C1=CC=CC=C1CC"
   obabel -:"$STYRENE" -ismi --gen2d -opng -O ${mollib}/inputs/STY.png
   obabel -:"$STYRENE" -ismi -omol2 -h --gen3d --title "Styrene-active" | \
      sed s/" 7 C "/" 7 C1"/ | \
      sed s/" 8 C "/" 8 C2"/ | \
      sed s/"UNL1"/"STY "/ > ${mollib}/inputs/STY.mol2

   # Components of BisGMA = HIE+GPA+HIE
   # Bisphenol-A
   PHENOL="C1=CC=C(O)C=C1"
   BPA="CC($PHENOL)($PHENOL)C"
   obabel -:"$BPA" -ismi --gen2d -opng -O ${mollib}/inputs/BPA.png -xp 600
   obabel -:"$BPA" -ismi --gen3d -h -omol2 --title "Bisphenol A" \
            | sed s/"UNL1"/"BPA "/ \
            | sed s/" 7 O "/" 7 O1"/ \
            | sed s/"14 O "/"14 O2"/ \
            > ${mollib}/inputs/BPA.mol2

   # 2-hydroxypropyl isopropyl ester
   HIE="CC(C(=O)OCC(O)C)C"
   obabel -:"$HIE" -ismi --gen2d -opng -O ${mollib}/inputs/HIE.png -xp 600
   obabel -:"$HIE" -ismi --gen3d -h -omol2 --title "2-hydroxypropyl isopropyl ester" \
            | sed s/"UNL1"/"HIE "/ \
            | sed s/"10 C "/"10 C2"/ \
            | sed s/" 2 C "/" 2 C1"/ \
            | sed s/" 7 C "/" 7 C3"/ \
            | sed s/" 9 C "/" 9 C4"/ \
            > ${mollib}/inputs/HIE.mol2

   htpolynet run -diag diagnostics.log GMASTY.yaml &> console.log

Note that we are actually generating three different monomeric species: styrene, bisphenol A (BPA), and 2-hydroxypropyl isopropyl ester (HIE).  This is because we are going to tell ``HTPolyNet`` to build the bisGMA molecule by reacting the esters onto the phenolic hydroxyls on BPA.

As in the other two tutorials, if you have a ``mol-env`` conda enviroment, just issue ``./run.sh`` to start the build.
