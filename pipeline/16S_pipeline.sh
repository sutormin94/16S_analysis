#!/bin/bash

python sorter.py

bash quality_control.sh

Rscript DADA2_part_1.R

bash mmseq_clustering.sh

### for removing temporary directories
#mkdir /home/niagara/Storage/MetaRus/V_Mamontov/16S_results/OTUs_output/WSBS_sponges

python mmseq2_clustering.py

Rscript DADA2_part_2.R

### for removing temporary directories
#rm -r /home/niagara/Storage/MetaRus/V_Mamontov/results_16S/mmseq
#rm -r /home/niagara/Storage/MetaRus/V_Mamontov/resutls_16S/Trimmed

echo "The job is completed!"

