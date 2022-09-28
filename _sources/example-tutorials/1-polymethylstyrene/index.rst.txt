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
   lib/  README.md  pMSTY.yaml  run.sh

The bash script ``run.sh`` is a *suggested* way to run the build, and it contains the commands to generate the input ``mol2`` file for activated methylstyrene (i.e., ethylmethylbenzene, or "EMB").  A listing of ``run.sh`` appears below:

.. code-block:: bash

   #!/bin/bash -l
   # A simple driver to demonstrate usage of HTPolyNet
   #
   # Cameron F. Abrams cfa22@drexel.edu

   # activate the conda environment we configured to run htpolynet
   conda activate mol-env

   sysname="pMSTY"
   mollib='./lib/molecules'

   # clean out any existing project directories and the library
   rm -rf proj-* *log ${mollib}/parameterized/*

   # Generate a 2d picture and a 3d mol2 file for methylstyrene's activated
   # form, ethylmethylbenzene (EMB) and rename reactive atoms
   echo "C1=CC(C)=CC=C1CC" | obabel -ismi --gen2d -opng -O ${mollib}/EMB.png
   echo "C1=CC(C)=CC=C1CC" | obabel -ismi --gen3d -h -omol2 --title "EMB" | \
         sed s/" 8 C "/" 8 C1"/ | \
         sed s/" 9 C "/" 9 C2"/ | \
         sed s/"UNL1"/"EMB "/ > ${mollib}/EMB.mol2

   # launch the build
   htpolynet run -diag diagnostics.log ${sysname}.yaml &> console.log

Note that we are activating the ``mol-env`` conda environment; that may not be necessary in your case, depending on how you configured your Python instance.  Note the two ``obabel`` commands: the first simply draws the chemical structure of the monomer as  PNG file, and the second generates the 3D structure in ``mol2`` format.  The ``sed`` commands are used to give unique atom names to atoms 8 and 9, and the default residue name ``UNL1`` is replaced with ``EMB``.  

So running the build just requires launching the script:

.. code-block:: bash

   $ ./run.sh &
   

The following sections guide the rest of the way through the tutorial.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   monomer
   reactions
   configuration
   run
   results
