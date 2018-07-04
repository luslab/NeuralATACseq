#!/bin/bash
#==============================================================================#
#     ATAC-seq signal in NMP and NMP/SC sites as region heatmap (Figure 4D)
#==============================================================================#
# Define Paths
BEDDIR=/home/steinhs/projects/atacseq_vicki/results/deseq2/mm10/
BWDIR=/home/steinhs/projects/atacseq_vicki/results/atac_conregions_KO/bw

# Compute summits from given region beds 
find -L $BEDDIR -type f | grep -E '(2_4|3_4)' |
xargs -I% sh -c 'NF=$(basename % .bed)"_summits.bed"; cat % |
awk "{mid=int((\$3-\$2)/2); print \$1,\$2+mid,\$2+mid+1}" | sed "s/ /\t/g" | bedtools sort -i - > $NF' 

# Compute matrix from *.bed and *.bw
BW=$(find $BWDIR | grep D3 | grep bw$ | sort)
computeMatrix reference-point -S $BW -R *.bed -a 2500 -b 2500 \
                              -out NMP_diffRegions.matrix.gz -p 15

# Plot heatmaps
LABEL=$(ls | grep bed$ | sort | cut -d"_" -f2)
SAMPLELABEL=$(find $BWDIR | grep D3 | grep bw$ | sort |
              xargs -I% basename % .bw | sed 's/_merged.*//g')

plotHeatmap  -m NMP_diffRegions.matrix.gz -out NMP_diffRegions_heatmap.pdf \
             --colorMap Reds --missingDataColor White --refPointLabel 'Summit' \
             --samplesLabel $SAMPLELABEL --regionsLabel $LABEL --xAxisLabel '' \
             --zMax 0.5 0.5 0.5 0.5

