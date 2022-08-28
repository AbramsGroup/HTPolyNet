.. _pms_tutorial:

Poly(4-methyl styrene)
======================

This tutorial illustrates how to build an all-atom model of a poly(methyl styrene) system using HTPolyNet.

The first step is to get to a clean working directory, and then use ``htpolynet fetch-example`` to setup this example.

.. code-block:: console

   $ mkdir my_pms
   $ cd my_pms
   $ htpolynet fetch-example -n 2
   $ ls
   2-polymethylstyrene/
   $ cd 2-polymethylstyrene
   $ ls
   lib/  pMSTY-hi.yaml  pMSTY-lo.yaml  run.sh

The bash script ``run.sh`` is a *suggested* way to run both a low-cure and high-cure build, and it contains the commands to generate the input ``mol2`` file for activated methylstyrene (i.e., ethylmethylbenzene, or "EMB").

.. code-block:: bash

   #!/bin/bash -l
   sysname="pMSTY"
   conda activate mol-env
   mollib='./lib/molecules/inputs'
   rm -rf proj-* *log lib/molecules/parameterized/*
   echo "C1=CC(C)=CC=C1CC" | obabel -ismi --gen2d -opng -O ${mollib}/EMB.png
   echo "C1=CC(C)=CC=C1CC" | obabel -ismi --gen3d -h -omol2 --title "EMB" | \
         sed s/" 8 C "/" 8 C1"/ | \
         sed s/" 9 C "/" 9 C2"/ | \
         sed s/"UNL1"/"EMB "/ > ${mollib}/EMB.mol2
   htpolynet run -diag diagnostics-lo.log ${sysname}-lo.yaml &> lo.log
   htpolynet run -diag diagnostics-hi.log ${sysname}-hi.yaml &> hi.log

So, if you have a ``mol-env`` conda enviroment, just issue ``./run.sh`` to start the builds.  The section "run" below walks through the build process, while the first three sections describe in detail how to create the monomer input and configuration file, with a special emphasis on explaining the reactions' syntax.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   monomer
   reactions
   configuration
   run
   results
