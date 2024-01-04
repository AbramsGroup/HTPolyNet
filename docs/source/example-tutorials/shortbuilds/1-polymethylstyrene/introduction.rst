.. _pms_introduction:

Introduction
------------

The first step is to get to a clean base directory, and then use ``htpolynet fetch-example`` to setup this example.

.. code-block:: console

   $ mkdir my_pms
   $ cd my_pms
   $ htpolynet fetch-example -n 2
   $ ls
   2-polymethylstyrene/
   $ tree 2-polymethylstyrene
   2-polymethylstyrene/
   ├── lib
   │   └── molecules
   │       ├── inputs
   │       └── parameterized
   ├── pMSTY.yaml
   ├── README.md
   └── run.sh

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
   
Results will appear in the files ``diagnostics.log`` and ``console.log``, and the project directory ``proj-0/``.

The first invocation of ``htpolynet run`` performs all needed parameterizations *and* runs the first system build.  These parameterizations appear in the library ``lib/molecules/parameterized/``.  All subsequent invocations of ``htpolynet run`` will use these parameterizations directly without re-performing them.

So, you could build multiple replicas by simply issuing more ``htpolynet run`` commands in series, for example:

.. code-block:: bash

   $ htpolynet run -diag diagnositics-1 pMSTY.yaml &> console-1.log
   $ htpolynet run -diag diagnositics-2 pMSTY.yaml &> console-2.log
   $ htpolynet run -diag diagnositics-3 pMSTY.yaml &> console-3.log

These will generate three new replicas in ``proj-1/``, ``proj-2/``, and ``proj-3/``, respectively.

Although this makes it seem quite simple to run ``htpolynet``, the aim of the tutorial is to show `how` it works to a level of detail sufficient to allow a user to develop their own new systems.  The remainder of this tutorial considers the major elements of ``htpolynet run`` and its usage.  Next we consider how to generate the :ref:`monomer inputs <tutorial_pms_monomer>`.