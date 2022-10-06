#!/bin/bash
source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate edta
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/GF1/11-edta/chr-edta/EDTA/bin/LTR_FINDER_parallel/bin
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/GF1/11-edta/chr-edta/EDTA/bin/LTR_HARVEST_parallel/bin
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/GF1/11-edta/chr-edta/EDTA/util
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/GF1/11-edta/chr-edta/EDTA
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/GF1/11-edta/chr-edta/EDTA/bin/HelitronScanner
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/GF1/11-edta/chr-edta/EDTA/bin/LTR_FINDER_parallel
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/GF1/11-edta/chr-edta/EDTA/bin/LTR_HARVEST_parallel
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/GF1/11-edta/chr-edta/EDTA/bin/README.md
PATH=$PATH:/work/LAS/mhufford-lab/arnstrm/GF1/11-edta/chr-edta/EDTA/bin/TIR-Learner2.5
cds=$2
genome=$1
EDTA.pl --genome $genome --species others --cds $cds --anno 1 --threads 36
