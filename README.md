# HTPolyNet
> High-Throughput Polymer Network Atomistic Simulations

HTPolyNet is a Python utility for generating atomistic models of cross-linked polymer networks together with appropriate topology and parameter files required for molecular dynamics simulations using Gromacs.  It is intended as a fully automated system builder requiring as inputs only the molecular structures of any monomer species, a description of the polymerization chemistry, and a handful of options describing desired system size and composition.  HTPolyNet uses the Generalized Amber Force Field for atom-typing and parameter generation.

## Installation

From PyPI:
```bash
pip install htpolynet
```

From source:
```bash
git clone git@github.com:AbramsGroup/HTPolyNet.git
cd HTPolyNet
pip install -e .
```

Once installed, the user has access to the main `htpolynet` command.

The programs ``antechamber``, ``parmchk2`` and ``tleap`` from AmberTools must be in your path.  These can be installed using the ``ambertools`` package from ``conda-forge`` or compiled from source.


## Release History

* 1.0.8
    * uses `chordless_cycles` to find rings
* 1.0.7.2
    * moved Library package to resources subpackage of HTPolyNet.HTPolyNet
* 1.0.6
    * gmx-style analyze subcommand added
* 1.0.5
    * Post-build MD simulations and plotting functionalities added
* 1.0.2
    * Enhanced molecule-network graph drawing in the plot subcommand
* 1.0.1
    * Fixed atom index assignment issue for systems with more than 100,000 atoms
* 1.0.0
    * First release
* 0.0.1
    * Initial beta version

## Meta

Cameron F. Abrams – cfa22@drexel.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

[https://github.com/AbramsGroup](https://github.com/AbramsGroup/)

## Contributing

1. Fork it (<https://github.com/AbramsGroup/HTPolyNet/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

