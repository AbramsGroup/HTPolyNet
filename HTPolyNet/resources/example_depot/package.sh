#!/bin/bash
DEPOT='/home/cfa/Git/HTPolyNet/Library/example_depot'
package () {
    tar --exclude="*/*/*/*/*" -zvcf ${DEPOT}/${1}.tgz ${1}/README.md ${1}/run.sh ${1}/*.yaml ${1}/lib/molecules/
}

# for r in 0-liquid-styrene \
#          1-polystyrene \
#          2-polymethylstyrene \
#          3-bisgma-styrene-thermoset \
#          4-pacm-dgeba-epoxy-thermoset \
#          5-dfda-fde-epoxy-thermoset \
#          6-htpb-ipdi; do 

#     package $r

# done
