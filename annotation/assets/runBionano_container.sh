#!/bin/bash
if [ "$#" -lt 2 ] ; then
echo "please provide:"
echo -e "\t\t(1) contigs/scaffolds in fasta format"
echo -e "\t\t(2) optical map in cmap format"
echo "";
echo "";
echo "contig/scaffold ids should be unique";
echo "cmap file is usually located in `output/contigs/exp_refineFinal1` folder";
echo "";
echo "./runBionano_container.sh.sh <contigs.fasta> <EXP_REFINEFINAL1.cmap>" ;
echo "";
exit 0;
fi
# load modules
ml purge
ml singularity
# set variables
image="/work/LAS/mhufford-lab/arnstrm/PanAnd/containers/bionano.v2.sif"
scmd="singularity exec --bind $PWD $image"
tdate=$(date +"%d%^b%Y")
genome=$1
cmap=$2
base="$(basename ${genome%.*})_$(basename ${cmap%.*})"
outdir=$(pwd)/$(basename ${genome%.*})_$(basename ${cmap%.*})
config="/opt/Solve3.7_10192021_74_1/HybridScaffold/1.0/hybridScaffold_DLE1_config.xml"
EBROOTBIONANOSOLVE="/opt/Solve3.7_10192021_74_1"
#output
mkdir -p ${outdir}
# run bionano solve
$scmd perl $EBROOTBIONANOSOLVE/HybridScaffold/1.0/hybridScaffold.pl \
-c ${config} \
-b ${cmap} \
-n ${genome} \
-u CTTAAG \
-z results_${base}_${tdate}.zip \
-w status_${base}_${tdate}.txt \
-B 2 \
-N 2 \
-g \
-f \
-r /opt/Solve3.7_10192021_74_1/RefAligner/1.0/sse/RefAligner \
-p /opt/Solve3.7_10192021_74_1/Pipeline/1.0 \
-o ${outdir}

