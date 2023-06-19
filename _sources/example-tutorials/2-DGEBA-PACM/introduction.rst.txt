.. _dgeba_pacm_introduction:

Introduction
============

The first step is to get to a clean working directory, and then use ``htpolynet fetch-example`` to setup this example.

.. code-block:: console

   $ mkdir my_dgeba_pacm
   $ cd my_dgeba_pacm
   $ htpolynet fetch-example -n 4
   $ ls
   4-pacm-dgeba-epoxy-thermoset/
   $ cd 4-pacm-dgeba-epoxy-thermoset
   $ ls
   DGE-PAC-hi.yaml  DGE-PAC-lo.yaml  lib/  run.sh

The bash script ``run.sh`` is a *suggested* way to run both a low-cure and high-cure build, and it contains the commands to generate the input ``mol2`` file for activated DGEBA and PACM.

.. code-block:: bash

   #!/bin/bash -l
   sysname="DGE-PAC"
   conda activate mol-env
   mollib='./lib/molecules/inputs'
   # DGEBA
   echo "CC(C)(C1=CC=C(C=C1)OCC(O)C)C3=CC=C(C=C3)OCC(O)C" | \
               obabel -ismi --gen2d -opng -O ${mollib}/DGE.png
   echo "CC(C)(C1=CC=C(C=C1)OCC(O)C)C3=CC=C(C=C3)OCC(O)C" | \
               obabel -ismi -h --gen3d -omol2 --title "DGE" | \
               sed s/"UNL1   "/"DGE    "/ | \
               sed s/"14 C "/"14 C1"/ | \
               sed s/"25 C "/"25 C2"/ | \
               sed s/"12 C "/"12 C3"/ | \
               sed s/"23 C "/"23 C4"/ | \
               sed s/"13 O "/"13 O1"/ | \
               sed s/"24 O "/"24 O2"/ > ${mollib}/DGE.mol2
   # PACM
   echo "C1CC(CCC1CC2CCC(CC2)N)N" | \
               obabel -ismi --gen2d -opng -O ${mollib}/PAC.png
   echo "C1CC(CCC1CC2CCC(CC2)N)N" | obabel -ismi -h --gen3d -omol2 --title "PAC" | \
               sed s/"UNL1   "/"PAC    "/ | \
               sed s/"14 N "/"14 N1"/ | \
               sed s/"15 N "/"15 N1"/ | \
               sed s/"3 C "/"3 C1"/ | \
               sed s/"11 C "/"11 C1"/ > ${mollib}/PAC.mol2
   rm -rf proj-* *log lib/molecules/parameterized/*
   htpolynet run -diag diagnostics-lo.log ${sysname}-lo.yaml &> lo.log
   htpolynet run -diag diagnostics-hi.log ${sysname}-hi.yaml &> hi.log

So, if you have a ``mol-env`` conda enviroment, just issue ``./run.sh`` to start the builds.  The section "run" below walks through the build process, while the first three sections describe in detail how to create the monomer input and configuration file, with a special emphasis on explaining the reactions' syntax.
