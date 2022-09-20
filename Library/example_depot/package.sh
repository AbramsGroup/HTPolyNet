#!/bin/bash
DEPOT='/home/cfa/Git/HTPolyNet/Library/example_depot'
for r in 0-liquid-styrene \
         1-polystyrene \
         2-polymethylstyrene \
         3-bisgma-styrene-thermoset \
         4-pacm-dgeba-epoxy-thermoset \
         5-dfda-fde-epoxy-thermoset \
         6-htpb-ipdi; do 

    tar --exclude="*/*/*/*/*" -zvcf ${DEPOT}/${r}.tgz ${r}/README.md ${r}/run.sh ${r}/*.yaml ${r}/lib/molecules/

done
