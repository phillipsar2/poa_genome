# Title when I come up with a paper name

Workflow for the phylogenetic species validation, analysis of genetic diversity, and evaluation of population strucutre of the *Poa pratensis* population panel. 

Analyses are partially implemented as a Snakemake workflow 
(https://snakemake.readthedocs.io/en/stable/index.html). 

There are three groups of analyses:
1. Generating consenus gene sequencing for *Poa* species confirmation
2. Population genetics of the *Poa* population panel
3. Determining ploidy of the *P. pratensis* genome

## Workflow

## Project Organization
 
├── README.md <- The top-level README.md for Poa enthusiasts that want to replicate these analyses.   
├── rules    
|   ├── consensus_seq.smk   
|   ├── pop_gen.smk   
|   └── ref_AB.smk   
├── environment.yml   
├── scripts   
│   ├── allelebalance_filter.sh   
│   └── allelebalance_filter.R   
├── notebooks   
│   └──    
├── data   
│   ├── raw <- The original WGS data dump.   
│   ├── interm  <- Intermediate data that has been transformed.   
│   ├── processed <- The final datasets for modeling.   
│   ├── vcf <- The unfiltered vcfs.   
│   ├── gene <- Genes downloaded from NCBI.    
│   └── genome <- The reference genome.   
├── reports <- Generated analyses as HTML, PDF, or .txt.    
├── Snakefile   
├── config.py   
├── results   
├── submit.json   
└── submit.sh   

