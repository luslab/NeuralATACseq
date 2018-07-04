#!/bin/bash
#==============================================================================#
#     Cdx2 ChIP-seq signal in SOM regions overlapping
#           Cdx2 peak as heatmap (Figure 5B)
#==============================================================================#
# Define paths
BEDDIR=/home/steinhs/projects/atacseq_vicki/results/cdx2/regions/mm10/

# Compute summits from given region beds 
find -L $BEDDIR -type f |grep bed | xargs -I% sh -c '
NF=$(basename % .bed)"_summits.bed";
cat % | awk "{mid=int((\$3-\$2)/2); print \$1,\$2+mid,\$2+mid+1}" |
sed "s/ /\t/g" | bedtools sort -i - > $NF' 

# Finde bigWigs for given histone TF
TF=Cdx2-FGF
BWDIR=/home/steinhs/projects/atacseq_vicki/data/chipseq/mazzoni2013_GSE39433/MNP
BWS=$(find $BWDIR -type f | grep bw$ | grep $TF | sort | grep subtr)
TFLABEL=$(echo $BWS | sed "s/ /\n/g" | xargs -I% basename % .bw |
          sed "s/_rep.*//g" | sed "s/EpiSCs//" | sed -E "s/^(-|_)//g")
echo $TFLABEL 

# Compute matrix from *.bed and *.bw
computeMatrix reference-point -S $BWS -R *.bed -a 5000 -b 5000\
                              -out $TF.matrix.gz -p 30

# Plot heatmaps
LABEL=$(find -L $BEDDIR -type f | sort | grep bed$ | xargs -I% basename % .bed | cut -d"_" -f1) 

plotHeatmap  -m $TF.matrix.gz -out $TF"_heatmap.png" --colorMap Reds \
             --missingDataColor White --refPointLabel 'Summit' \
             --samplesLabel $TFLABEL --regionsLabel $LABEL --xAxisLabel '' --zMin 0 
