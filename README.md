# HTPolyNet
Development of High-Throughput Polymer Network Atomistic Simulation

# Required packages
* [AmberTools19](https://ambermd.org/GetAmber.php#ambertools)
  - Options:
     * conda: `conda install -c ambermd ambertools` installs precompiled executables for ambertools into your Anaconda environment
     * compile from source:
       ```
       tar jxf AmberTools20.tar.bz2
       cd amber_src
       ./configure gnu
       source amber.sh
       make install
       ```
  - The executables needed specfically are `antechamber`, `parmchk2`, and `tleap`

* gromacs v2016 (Don't choose v2018)
  - Notes for 2016.6 on OpenSUSE Leap 15.2 with cuda 11 and openMPI:
    ```
    tar xfz gromacs-2016.6.tar.gz
    cd gromacs-2016.6
    mkdir build
    cd build
    cmake ..  -DGMX_BUILD_OWN_FFTW=ON -DREGRESSIONTEST_DOWNLOAD=ON -DGMX_MPI=on -DGMX_GPU=on -DCMAKE_INSTALL_PREFIX=/usr/local/gromacs-2016.6
    make
    make check
    sudo make install
    source /usr/local/gromacs-2016.6/bin/GMXRC  (added to ~/.bashrc)
    ```
  

* python3.x (I'm using Anaconda)
  - numpy
  - pandas
  - scipy
  - parmed
    ```
    conda install -c omnia parmed
    ```
  - GromacsWrapper (there is no conda package)
    ```
    pip install GromacsWrapper
    ````
    - Requires a custom `~/.gromacswrapper.cfg` file; to create a template:
      ```python
      >>> import gromacs
      >>> gromacs.config.setup()
      ```
    - Edit `~/.gromacswrapper.cfg`:
      ```
      release = 2016.6 # Type your GMX VERSION
      gmxrc = /usr/local/gromacs/bin/GMXRC # Type where you gmx execute file's location
      tools = gmx gmx_mpi
      ```