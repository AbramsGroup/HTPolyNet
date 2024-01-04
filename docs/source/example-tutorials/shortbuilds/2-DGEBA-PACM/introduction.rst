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
   $ tree 4-pacm-dgeba-epoxy-thermoset
   4-pacm-dgeba-epoxy-thermoset/
   ├── DGEPAC.yaml
   ├── lib
   │   └── molecules
   │       ├── inputs
   │       └── parameterized
   ├── README.md
   └── run.sh


The bash script ``run.sh`` is a *suggested* way to run a build, and it contains the commands to generate the input ``mol2`` files for activated DGEBA and PACM.

.. code-block:: bash

   #!/bin/bash -l
   # A simple driver to demonstrate usage of HTPolyNet
   #
   # Cameron F. Abrams cfa22@drexel.edu

   # activate the conda environment we configured to run htpolynet
   conda activate mol-env

   sysname="DGEPAC"
   mollib='./lib/molecules'

   # clean out any existing project directories and the library
   rm -rf proj-* *log ${mollib}/parameterized/*

   # Monomer structures:
   # DGEBA
   DGEBA="CC(C)(C1=CC=C(C=C1)OCC(O)C)C3=CC=C(C=C3)OCC(O)C"
   obabel -:"$DGEBA" -ismi --gen2d -opng -O ${mollib}/inputs/DGE.png
   obabel -:"$DGEBA" -ismi -h --gen3d -omol2 --title "DGE" | \
               sed s/"UNL1   "/"DGE    "/ | \
               sed s/"14 C "/"14 C1"/ | \
               sed s/"25 C "/"25 C2"/ | \
               sed s/"12 C "/"12 C3"/ | \
               sed s/"23 C "/"23 C4"/ | \
               sed s/"13 O "/"13 O1"/ | \
               sed s/"24 O "/"24 O2"/ > ${mollib}/inputs/DGE.mol2
   # PACM
   PACM="C1CC(CCC1CC2CCC(CC2)N)N" 
   obabel -:"$PACM" -ismi --gen2d -opng -O ${mollib}/inputs/PAC.png
   obabel -:"$PACM" -ismi -h --gen3d -omol2 --title "PAC" | \
               sed s/"UNL1   "/"PAC    "/ | \
               sed s/"14 N "/"14 N1"/ | \
               sed s/"15 N "/"15 N1"/ | \
               sed s/"3 C "/"3 C1"/ | \
               sed s/"11 C "/"11 C1"/ > ${mollib}/inputs/PAC.mol2

   # launch the build
   htpolynet run -diag diagnostics.log ${sysname}.yaml &> console.log

So, if you have a ``mol-env`` conda enviroment, just issue ``./run.sh`` to start the builds.  The section "run" below walks through the build process, while the first three sections describe in detail how to create the :ref:`monomer inputs <dgeba_pacm_monomers>` and :ref:`configuration file <dgeba_configuration_file>`, with a special emphasis on explaining the :ref:`syntax of the reactions <dgeba_reaction_dictionaries>`.
