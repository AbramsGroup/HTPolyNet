# HTPolyNet
Development of High-Throughput Polymer Network Atomistic Simulation

# Needed package
AmberTools
 - antechamber
 - parmck2
 - tleap

gromacs v2016 (Don't choose v2018)
python3
 - numpy
 - pandas
 - scipy
 - parmed
 - GromacsWrapper
   - To use this gromacs module, need to configure the module first after installation
   - Under the python3 console, import the gromacs module and use "gromacs.config.setup()" to generate the configure files (.gromacswrapper.cfg, .gromacswrapper)
   - modified file .gromacswrapper.cfg: # basic option
     - release = 2016.6 # Type your GMX VERSION
     - gmxrc = /usr/local/gromacs/bin/GMXRC # Type where you gmx execute file's location
     - tools = gmx gmx_mpi
