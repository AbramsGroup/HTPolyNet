Installation and Prerequisites
==============================

Installation
------------

Because HTPolyNet is under active development, it is best to clone the repository from Github and then 
pip install as in:

.. code-block:: console

    $ git clone git@github.com:AbramsGroup/HTPolyNet.git
    $ cd HTPolyNet
    $ pip install -e .

Software Prequisites
--------------------

Because HTPolyNet uses the Generalized Amber Force Field (GAFF) and Gromacs, these should both be installed, and the following executables should be in your path:

1. ``antechamber``
2. ``parmchk2``
3. ``tleap``
4. ``gmx`` or ``gmx_mpi``

Suggested approaches for obtaining these are below.

* `AmberTools <https://ambermd.org/GetAmber.php#ambertools>`_. You have at least two options:

    * conda:  ``conda install -c conda-force ambertools`` installs precompiled executables into your conda environment
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
   * compile from source:

     .. code-block:: console

         $ tar xfz gromacs-2022.1.tar.gz
         $ cd gromacs-2022.1
         $ mkdir build
         $ cd build
         $ cmake ..  -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_MPI=on -DGMX_GPU=on -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs-2022.1
         $ make
         $ make check
         $ sudo make install

     And add to ``~/.bashrc``:

     .. code-block:: console

         source /usr/local/gromacs-2022.1/bin/GMXRC
     
     This should provide access to the ``gmx`` command.  If you compiled an MPI version, you will instead generated ``gmx_mpi``; either of these commands can be used by HTPolyNet.



