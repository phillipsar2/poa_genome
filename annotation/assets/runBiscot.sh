#!/bin/bash

if [ "$#" -lt 1 ] ; then
echo "";
echo "";
echo -n "please provide: "
echo -e "genome contigs in fasta format"
echo "";
echo "runBiscot.sh  <genome>"
echo "";
exit 0;
fi


source /work/LAS/mhufford-lab/arnstrm/miniconda/etc/profile.d/conda.sh
conda activate biscot

cmapQ=$(find $(pwd) -name "EXP_REFINEFINAL1_bppAdjust_cmap_*_NGScontigs_HYBRID_SCAFFOLD_q.cmap")
if [ -z ${cmapQ+x} ]; then echo "query cmap file not found"; exit 1; fi

cmapR=$(find $(pwd) -name "EXP_REFINEFINAL1_bppAdjust_cmap_*_NGScontigs_HYBRID_SCAFFOLD_r.cmap")
if [ -z ${cmapR+x} ]; then echo "reference cmap file not found"; exit 1; fi

xmap=$(find $(pwd) -name "EXP_REFINEFINAL1_bppAdjust_cmap_*_NGScontigs_HYBRID_SCAFFOLD.xmap")
if [ -z ${xmap+x} ]; then echo "xmap file not found"; exit 1; fi

key=$(find $(pwd) -name "*_CTTAAG_0kb_0labels_key.txt.cut.txt")
if [ -z ${key+x} ]; then echo "key file not found"; exit 1; fi


genome=$1

biscot --cmap-ref $cmapR \
--cmap-1 $cmapQ \
--xmap-1 $xmap \
--key $key \
--contigs ${genome} \
--output biscot

