#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]; then
  echo "Usage: ./mol2lammps.sh Name.mol2 [charge] [integer/skip]";
else
  echo "Execute this script! Don't source!"; return
fi

DEBUG=1 #1 will retain intermediate files

#Usage Eg 1: source mol2lammps.sh Name.mol2 
#Usage Eg 2: source mol2lammps.sh Name.mol2 charge -1
#Usage Eg 3: source mol2lammps.sh Name.mol2 charge skip

Main ()
{
  CheckArguments $1 $2 $3
  CheckExecutables 
  Mol2toLAMMPS $1 $2 $3
  GenerateInputScript $1
}

CheckArguments() {
  if [ -z "$1" ]; then
    echo "Error: Please provide the name of mol2 file after script name!" & exit 1
  elif [ ! -f $1 ]; then
    echo "Error: mol2 file does not exist in current directory" & exit 1
  fi

  if [ ! -z "$2" ] && [ "$2" != "charge" ]; then
      echo "Error: Second keyword is optional, but must be charge" & exit 1; fi

  if [ "$2" = "charge" ]; then 
    if [ -z "$3" ]; then echo "Error: charge keyword must be followed by integer or keyword skip" & exit 1; fi
    if [[ "$3" =~ ^-?[0-9]+$ ]]; then echo "Info: Net charge specified!";
    elif [ "$3" = "skip" ]; then echo "Info: Charges will not be calculated";  
    else echo "Error: Third argument must be integer or skip" & exit 1; fi
  fi
}

CheckExecutables() {
  if [ ! $(command -v antechamber) &> /dev/null ]; then echo "Error: antechamber not found" & exit 1; fi
  if [ ! $(command -v parmchk2) &> /dev/null ]; then echo "Error: parmchk2 not found" &  exit 1; fi
  if [ ! $(command -v tleap) &> /dev/null ]; then echo "Error: tleap not found" & exit 1; fi
  if [ ! $(command -v python3) &> /dev/null ]; then echo "Error: python3 not found" & exit 1; fi
  if [ ! -f amber2lammps.py3 ]; then echo "Error: amber2lammps.py3 does not exist in current directory" & exit 1; fi
}

Mol2toLAMMPS () {
  BaseName=${1%.mol2}

  if [ -z "$2" ]; then antechamber -i $1 -fi mol2 -o ${BaseName}bcc.mol2 -fo mol2 -c bcc -pf y
  elif [ "$3" = "skip" ]; then antechamber -i $1 -fi mol2 -o ${BaseName}bcc.mol2 -fo mol2 -c dc -pf y
  else antechamber -i $1 -fi mol2 -o ${BaseName}bcc.mol2 -fo mol2 -c bcc -pf y -nc $3; fi
  
  parmchk2 -i ${BaseName}bcc.mol2 -f mol2 -o ${BaseName}.frcmod -a Y

  printf "loadamberparams %s.frcmod\n" ${BaseName} > leaprc.${BaseName}
  printf "a=loadmol2 %sbcc.mol2\n" ${BaseName} >> leaprc.${BaseName}
  printf "saveamberparm a %s.top %s.crd\n" ${BaseName} ${BaseName} >> leaprc.${BaseName}
  printf "quit" >> leaprc.${BaseName}
  tleap -f leaprc.${BaseName}
  python3 amber2lammps.py3 ${BaseName}

  if [ "${DEBUG}" -eq 0 ]; then
    rm ${BaseName}.top ${BaseName}.crd leaprc.${BaseName} sqm.in sqm.out sqm.pdb ${BaseName}bcc.mol2 ANTECHAMBER.FRCMOD
  fi
}

GenerateInputScript() {
  BaseName=${1%.mol2}
  cat > in.${BaseName} <<EndInput
units  real
neigh_modify  delay 5 every 1
atom_style  full
pair_style  lj/cut/coul/long 9.0
pair_modify  mix arithmetic tail yes
bond_style  harmonic
angle_style  harmonic
dihedral_style  charmm
improper_style  cvff
kspace_style  pppm 1e-04
kspace_modify  diff ad
boundary  p p p
special_bonds  amber
  
read_data Data.${BaseName}
  
minimize 0 0 10000 10000
  
replicate ${NX} ${NY} ${NZ}
  
variable MassG equal mass(all)/6.02214e23
variable Volcm3 equal vol*1.0e-24
variable Density equal v_MassG/v_Volcm3
  
thermo_style  custom step temp press vol etotal epair emol evdwl ecoul elong ebond eangle enthalpy pe ke v_Density lx ly lz cpu
thermo  1000

fix  MOMENTUM all momentum 100000 linear 1 1 1 angular
fix  SHAKE all shake 0.00001 25 0 b 4 7 8 10 11 15 17 a 3 6 10 13 14 15 16 17 18 19 25 27 28 32 m 1.008
fix  NPT all npt temp ${Temperature} ${Temperature} 100 iso 1.0 1.0 1000
fix  RECENTER all recenter INIT INIT INIT units box
  
run 10000
EndInput
}

Main $1 $2 $3
