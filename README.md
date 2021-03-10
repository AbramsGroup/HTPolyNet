# HTPolyNet
Development of High-Throughput Polymer Network Atomistic Simulation

# Required packages
* [AmberTools20](https://ambermd.org/GetAmber.php#ambertools)
  - antechamber
  - parmck2
  - tleap

* gromacs v2016 (Don't choose v2018)
  - Note, if you have cuda 11 or higher you may need to comment out the line that enables support for `compute_30` in `cmake/gmxManageNvccConfig.cmake`:
    ```
    #list (APPEND GMX_CUDA_NVCC_GENCODE_FLAGS "-gencode;arch=compute_30,code=sm_30")
    ```

* python3.x (I'm using Anaconda)
  - numpy
  - pandas
  - scipy
  - parmed  (`conda install -c omnia parmed`)
  - GromacsWrapper (`pip install GromacsWrapper`; there is no conda package)
    - To use this gromacs module, need to configure the module first after installation
    - Under the python3 console, import the gromacs module and use "gromacs.config.setup()" to generate the configure files (.gromacswrapper.cfg, .gromacswrapper)
    -  modified file .gromacswrapper.cfg: # basic option
      - release = 2016.6 # Type your GMX VERSION
      - gmxrc = /usr/local/gromacs/bin/GMXRC # Type where you gmx execute file's location
      - tools = gmx gmx_mpi
