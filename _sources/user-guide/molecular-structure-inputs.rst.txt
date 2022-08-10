.. _molecular_structure_inputs:

Molecular Structure Inputs
--------------------------

The intention of HTPolyNet is to generate amorphous, polymerized all-atom system models when provided 3D structure/topology information about any monomers and a description of the polymerization chemistry.  Any molecule you wish to use as a **monomer** must have a either a ``pdb``-format file or a ``mol2``-format file.  The basename of the file should be the name of the molecule.  For example, ``STY.mol2`` might be a styrene monomer.  This basename/monomer name is important; it is how the monomer is called forever inside HTPolyNet.  Input ``mol2`` and ``pdb`` files are all processed via antechamber as part of the parameterization algorithm.

Sample mol2 files for a few monomers are provided in the Library subpackage: ``Library/molecules/sample-inputs``.  You will likely want to create your own.  Most chemical structure drawing programs will output mol2 files.  OpenBabel can also generate them from SMILES strings.  The mol2 format HTPolyNet requires for a monomer is minimal and requires only three sections: ``@<TRIPOS>MOLECULE``, ``@<TRIPOS>ATOM`` and ``@<TRIPOS>BOND``.  No matter how you generate a mol2 file, if it corresponds to a molecule that can react, you must edit the mol2 file to give **unique atom names** to any reactive atoms or chiral carbons.  You will then use those unique names to refer specifically to those atoms.  The best way to identify these atoms is by inspecting the 3d structure of the molecule.

