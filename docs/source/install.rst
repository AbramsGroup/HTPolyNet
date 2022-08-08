Installation and Prerequisites
==============================

Python
------

If you use conda/anaconda, we recommended that you create a separate python environment running HTPolyNet:

.. code-block:: console

    $ conda create --name mol-env python
    $ conda activate mol-env


Once this environment is created and activated, you can install some prequisites:

.. code-block:: console

    $ conda install numpy scipy pandas ambertools


Installation
------------

Because HTPolyNet is under active development, we recommend installing from a freshly cloned Github repository:

.. code-block:: console

    $ git clone git@github.com:AbramsGroup/HTPolyNet.git
    $ cd HTPolyNet
    $ pip install -e .

If you created the recommended python environment, make sure it is activated before running ``pip install``!

Software Prequisites
--------------------

Because HTPolyNet uses the Generalized Amber Force Field (GAFF) and Gromacs, it expects the following executables are in your path:

1. ``antechamber``, ``parmchk2``, and ``tleap`` (from AmberTools)
2. ``gmx`` or ``gmx_mpi``

Suggested approaches for obtaining these are below.

* `AmberTools <https://ambermd.org/GetAmber.php#ambertools>`_. You have at least two options:

    * conda:  ``conda install ambertools`` installs precompiled executables into your active environment (note that if you created the separate environment as described above, you have already installed ambertools)
    * compile from source (requires ``csh``, ``flex``, and ``bison``):

    .. code-block:: console

        $ tar jxf AmberTools21.tar.bz2
        $ cd amber_src
        $ ./configure --no-X11 --skip-python gnu
        $ source amber.sh
        $ make install

    Either of these approaches should provide access to ``antechamber``, ``parmchk2``, and ``tleap``.

* `Gromacs <https://manual.gromacs.org/documentation/current/index.html>`_.  You likely have at least two options:

   * your linux package manager
   * compile from source.  This is a good option if you want to use GPUs and/or MPI.  The example below builds Gromacs with CUDA but without MPI (assuming you have CUDA installed):

     .. code-block:: console

         $ tar xfz gromacs-2022.1.tar.gz
         $ cd gromacs-2022.1
         $ mkdir build
         $ cd build
         $ cmake ..  -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_GPU=CUDA -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs
         $ make
         $ make check
         $ sudo make install

     And add to your ``~/.bashrc``:

     .. code-block:: console

         source /usr/local/gromacs/bin/GMXRC
     
     This should provide access to the ``gmx`` command.  If you additionally compiled an MPI version (using ``-DGMX_MPI=on`` in the ``cmake`` command), you will also have access to ``gmx_mpi``; either of these commands can be used by HTPolyNet.  Note that Gromacs 2016 and below have a version of ``gmx distance`` that limits the number of distances that can be calculated, so we (always) recommend using the latest Gromacs.

Finally, in the tutorials provided, we demonstrate the use of ``obabel`` from `OpenBabel <https://openbabel.org/wiki/Main_Page>`_ to generate initial molecular structure files in Sybyl MOL2 format from SMILES strings; though it is not strictly necessary, it is fairly convenient to use for this purpose.