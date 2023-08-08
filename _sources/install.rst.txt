##############################
Installation and Prerequisites
##############################

Software Prequisites
--------------------

The following commands should be in your path:

1. ``antechamber``, ``parmchk2``, and ``tleap`` (`AmberTools <https://ambermd.org/GetAmber.php#ambertools>`_, version 22 or higher); preferred installation via ``conda``.
2. ``gmx`` or ``gmx_mpi`` (`Gromacs <https://manual.gromacs.org/documentation/current/index.html>`_, version 2022.1 or higher); preferred installation via compiling from source.
3. ``obabel`` (`OpenBabel <https://github.com/openbabel/openbabel>`_); preferred installation via Linux distribution package.

If you use conda/anaconda, we recommended that you create a separate Python environment running ``HTPolyNet``:

.. code-block:: console

    $ conda create --name mol-env python
    $ conda activate mol-env

Once this environment is created and activated, you can optionally install some prequisites:

.. code-block:: console

    $ conda install numpy scipy pandas ambertools

(``ambertools`` is available from the ``conda-forge`` channel.)  If you don't do this, ``setup.cfg`` in ``HTPolyNet`` will do it for you, which might take a bit.

Installation
------------

The current stable version of ``HTPolyNet`` is available at PyPI:

.. code-block:: console

    $ pip install htpolynet

To install a development version of ``HTPolyNet`` you can instal from a freshly cloned Github repository:

.. code-block:: console

    $ git clone git@github.com:AbramsGroup/HTPolyNet.git
    $ cd HTPolyNet
    $ pip install -e .

If you created the recommended python environment, make sure it is activated before running ``pip install``!


Notes
-----

If you prefer to use more recent versions of AmberTools, Gromacs, or OpenBabel than your system currently provides, you can compile the latest versions from source.  **It is recommended that you deactivate any conda environment before performing any of these compilations.**

Compilation of AmberTools
#########################

Compilation of AmberTools requires ``csh``, ``flex``, and ``bison``:

.. code-block:: console

    $ tar jxf AmberTools21.tar.bz2
    $ cd amber_src
    $ ./configure --no-X11 --skip-python gnu
    $ source amber.sh
    $ make install

Compilation of Gromacs
######################

You can also compile Gromacs from source, if your Linux distibution doesn't include it in its package management, or you are on a big supercomputer.  The example below builds Gromacs with CUDA but without MPI (assuming you have CUDA installed):

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

Compilation of ``obabel``
#########################

If your system does not have ``obabel`` installed and your Linux distribution doesn't offer a package for it (or you are not root!), you can compile it from source.  Be sure to unpack `Eigen <https://eigen.tuxfamily.org/index.php?title=Main_Page>`_ first so that the ``conformer`` plug-in for ``obabel`` will work.  Below I demonstrate a session in which both the Eigen and OpenBabel source packages are downloaded to ``~/Downloads`` and are unpacked in the directory ``~/build/``, and the OpenBabel installation directory is ``~/opt/obabel``.

.. code-block:: console

    $ cd ~/build
    $ tar jxf ~/Downloads/eigen-3.4.0.tar.bz2  # unpack only -- no need to compile
    $ tar jxf ~/Downloads/openbabel-3.1.1.tar.bz2
    $ cd openbabel-3.1.1
    $ mkdir build
    $ cd build
    $ cmake .. -DEIGEN3_INCLUDE_DIR=${HOME}/build/eigen-3.4.0/ -DCMAKE_INSTALL_PREFIX=${HOME}/opt/obabel
    $ make
    $ make test
    $ make install

You will need to ensure that ``${HOME}/opt/babel/bin`` is in your ``PATH``, ``${HOME}/opt/babel/lib`` is in your ``LD_LIBRARY_PATH``, and that the environment variable ``BABEL_LIBDIR`` is set to ``${HOME}/opt/babel/lib``.

Other Prequisites
-----------------

In order to use ``HTPolyNet`` effectively, it is recommended that you have good working knowledge of the following:

1. MD simulation in general and Gromacs specifically;
2. the General Amber Force Field (GAFF), including in particular

   a. how to use ``antechamber``, ``tleap``, and ``parmchk2`` to generate GAFF parameterizations; and
   b. how to use these parameterizations inside Gromacs; and

3. Polymer chemistry, at least for the systems you are interested in simulating.
