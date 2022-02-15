Workflow for the phylogenetic species validation, analysis of genetic diversity, and evaluation of population structure of the *Poa pratensis* population panel. 

Analyses are partially implemented as a Snakemake workflow 
(https://snakemake.readthedocs.io/en/stable/index.html). 

There are three groups of analyses broken into seperate rule files:
1. Generating consenus gene sequencing for *Poa* species confirmation (rules/consensus_seq.smk)
2. Population genetics of the *Poa* population panel (rules/pop_gen.smk)
3. Determining ploidy of the *P. pratensis* genome (rules/ref_AB.smk & scripts/predict_chromcount.R)

## Project Organization
<pre>
├── README.md <- The top-level README.md for Poa enthusiasts that want to replicate these analyses.   
├── rules    
|   ├── consensus_seq.smk   
|   ├── pop_gen.smk   
|   └── ref_AB.smk   
├── environment.yml   
├── scripts   
│   ├── allelebalance_filter.sh   
│   └── allelebalance_filter.R   
├── data   
│   ├── raw <- The original WGS data dump.   
│   ├── interm  <- Intermediate data that has been transformed.   
│   ├── processed <- The final datasets.
│   ├── vcf <- The unfiltered vcfs.   
│   ├── gene <- Genes downloaded from NCBI for phylogenetic analyses.    
│   └── genome <- The reference genome.   
├── reports <- Generated analyses as HTML, PDF, or .txt.    
├── Snakefile   
├── config.py   
├── submit.json   
└── submit.sh   
</pre>

