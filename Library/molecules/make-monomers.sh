#!/bin/bash
#
# HTPolyNet
#
# Examples of using obabel and SMILES to generate monomer structures
# ready to use
#
# Cameron F. Abrams, cfa22@drexel.edu

# STYRENE (ethylstyrene)
echo "C1=CC=CC=C1CC" | obabel -ismi --gen2d -opng -O pics/STY.png -xp 600 -xt
echo "C1=CC=CC=C1CC" | obabel -ismi -omol2 -h --gen3d --title STY | \
     sed s/" 7 C "/" 7 C1"/ | \
     sed s/" 8 C "/" 8 C2"/ | \
     sed s/"UNL1   "/"STY    "/ > inputs/STY.mol2

# METHYLSTYRENE (ethylmethylbenzene)
echo "C1=CC(C)=CC=C1CC" | obabel -ismi --gen2d -opng -O pics/EMB.png -xp 600 -xt
echo "C1=CC(C)=CC=C1CC" | obabel -ismi --gen3d -h -omol2 --title "EMB" | \
      sed s/" 8 C "/" 8 C1"/ | \
      sed s/" 9 C "/" 9 C2"/ | \
      sed s/"UNL1"/"EMB "/ > inputs/EMB.mol2

# BisGMA
PHENOLTHING="C1=CC=C(OCC(O)COC(=O)C(C)C)C=C1"
obabel -:"C(C)(C)($PHENOLTHING)($PHENOLTHING)" -ismi -opng -O pics/GMA.png -xp 600 -xt
obabel -:"C(C)(C)($PHENOLTHING)($PHENOLTHING)" -ismi --gen3d -h -omol2 --title GMA | \
      sed s/"16 C "/"16 C1"/ | \
      sed s/"33 C "/"33 C2"/ | \
      sed s/"17 C "/"17 C3"/ | \
      sed s/"34 C "/"33 C4"/ | \
      sed s/"10 C "/"10 C5"/ | \
      sed s/"27 C "/"27 C6"/ | \
      sed s/"UNL1   "/"GMA    "/ > inputs/GMA.mol2

# DGEBA
obabel -:"CC(C)(C1=CC=C(C=C1)OCC(O)C)C3=CC=C(C=C3)OCC(O)C" -ismi --gen2d -opng -O pics/DGE.png -xp 600 -xt
obabel -:"CC(C)(C1=CC=C(C=C1)OCC(O)C)C3=CC=C(C=C3)OCC(O)C" -ismi -h --gen3d -omol2 --title "DGE" | \
            sed s/"UNL1   "/"DGE    "/ | \
            sed s/"14 C "/"14 C1"/ | \
            sed s/"25 C "/"25 C2"/ | \
            sed s/"12 C "/"12 C3"/ | \
            sed s/"23 C "/"23 C4"/ | \
            sed s/"13 O "/"13 O1"/ | \
            sed s/"24 O "/"24 O2"/ > inputs/DGE.mol2

# PACM
obabel -:"C1CC(CCC1CC2CCC(CC2)N)N" -ismi --gen2d -opng -O pics/PAC.png -xp 600 -xt
obabel -:"C1CC(CCC1CC2CCC(CC2)N)N" -ismi -h --gen3d -omol2 --title "PAC" | \
            sed s/"UNL1   "/"PAC    "/ | \
            sed s/"14 N "/"14 N1"/ | \
            sed s/"15 N "/"15 N1"/ | \
            sed s/"3 C "/"3 C1"/ | \
            sed s/"11 C "/"11 C1"/ > inputs/PAC.mol2

# FDE and DFDA
FOURAMINOFURANYL="C1OC(CN)=CC=1"
FURANYL=c1occc1
DFDA_SMILES="C(${FOURAMINOFURANYL})${FOURAMINOFURANYL}"
FDE_SMILES="N((C${FURANYL})CC(O)C)CC(O)C"

# DFDA
echo $DFDA_SMILES | obabel -ismi --gen2d -opng -O pics/DFA.png -xp 600 -xt
echo $DFDA_SMILES | obabel -ismi -h --gen3d --title DFA -opdb |
sed s/"UNL "/"DFA "/ | \
sed s/"HETATM    6  N "/"HETATM    6  N1"/ | \
sed s/"HETATM   13  N "/"HETATM   13  N2"/ > inputs/DFA.pdb

# FDE
echo $FDE_SMILES | obabel -ismi --gen2d -opng -O pics/FDE.png -xp 600 -xt
echo $FDE_SMILES | obabel -ismi -h --gen3d --title FDE -opdb | \
sed s/"UNL "/"FDE "/ |
sed s/"HETATM   11  C "/"HETATM   11  C1/" |
sed s/"HETATM   15  C "/"HETATM   15  C2/" |
sed s/"HETATM    9  C "/"HETATM    9  C3/" |
sed s/"HETATM   13  C "/"HETATM   13  C4/" > inputs/FDE.pdb