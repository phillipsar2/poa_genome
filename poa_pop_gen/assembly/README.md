# Genome assembly

## Long reads error correction

Falcon was used for correcting the PacBio Subreads. The subreads in bam format was first converted to fasta format and config file was created to run Falcon.

```bash
for bam in *.bam; do
  samtools fasta --threads 36 $bam > ${bam%.*}.fasta;
done
```
The config file for Falcon is located [here](assets/falcon.cfg)

Falcon was kept running on the cluster using the script [`runFalconEC.sh`](runFalconEC.sh)

Upon completion, the error-corrected reads were concatenated to a single file as follows:

```bash
for cns in 0-rawreads/cns-runs/cns_00???/uow-00/consensus.cns_00???.fasta; do
  cat $cns >> consensus.merged.fa;
done
gzip consensus.merged.fa
```

The final report for the falcon is as follows:

```json
{
    "genome_length": 6400000000,
    "length_cutoff": 19801,
    "preassembled_bases": 104518236733,
    "preassembled_coverage": 16.331,
    "preassembled_esize": 21827.09,
    "preassembled_mean": 16985.699,
    "preassembled_n50": 20457,
    "preassembled_p95": 33979,
    "preassembled_reads": 6153308,
    "preassembled_seed_fragmentation": -1.0,
    "preassembled_seed_truncation": -1.0,
    "preassembled_yield": 0.408,
    "raw_bases": 397639775813,
    "raw_coverage": 62.131,
    "raw_esize": 28217.229,
    "raw_mean": 16309.184,
    "raw_n50": 25617,
    "raw_p95": 43980,
    "raw_reads": 24381341,
    "seed_bases": 256008216776,
    "seed_coverage": 40.001,
    "seed_esize": 37180.216,
    "seed_mean": 32960.457,
    "seed_n50": 33925,
    "seed_p95": 56114,
    "seed_reads": 7767132
}
```

Version: `falcon-kit 1.8.1`, `pypeflow 2.3.0`


## Assembly

Canu was used for assembling the error-corrected PacBio Reads. The input file `consensus.merged.fa.gz` was run with the script [`runCanu.slurm`](runCanu.slurm), using the config file [`canu.cfg`](assets/canu.cfg).

Contig stats for Canu is as follows:

```json
{
    "Assumed genome size (Mbp)": 3330.00,
    "Number of scaffolds": 27953,
    "Total size of scaffolds": 6243783987,
    "Total scaffold length as percentage of assumed genome size": 187.5,
    "Longest scaffold": 6812402,
    "Shortest scaffold": 1009,
    "Number of scaffolds > 1K nt": 27953,
    "Number of scaffolds > 10K nt": 24818,
    "Number of scaffolds > 100K nt": 11916,
    "Number of scaffolds > 1M nt": 1256,
    "Number of scaffolds > 10M nt": 0,
    "Mean scaffold size": 223367,
    "Median scaffold size": 65519,
    "N50 scaffold length": 642088,
    "L50 scaffold count": 2732,
    "NG50 scaffold length": 1118863,
    "LG50 scaffold count": 986,
    "N50 scaffold - NG50 scaffold length difference": 476775,
    "scaffold %A": 27.25,
    "scaffold %C": 22.74,
    "scaffold %G": 22.75,
    "scaffold %T": 27.26,
    "scaffold %N": 0.00,
    "scaffold %non-ACGTN": 0.00,
    "Number of scaffold non-ACGTN nt": 0
}
```

Version: `Canu 2.0`

## Hybrid scaffolding

BioNano Solve was used for hybrid scaffolding. The scaffolding was done using the [`runBionano.sh`](runBionano.sh) script. The summary stats are as follows:


| FileType                                              | Count   | Min length (Mbp)  | Median length (Mbp)  | Mean length (Mbp)  | N50 length (Mbp)  | Max length (Mbp)  | Total length (Mbp)  |
|-------------------------------------------------------|--------:|------------------:|---------------------:|-------------------:|------------------:|------------------:|--------------------:|
| Original BioNano Genome Map               |     183 |             0.078 |               30.347 |             34.259 |            61.972 |           178.078 |           6,269.426 |
| Bpp-adjusted BioNano Genome Map         |     183 |             0.077 |               29.897 |             33.751 |            61.052 |           175.434 |           6,176.352 |
| Original NGS sequences                    |  27,953 |             0.001 |                0.066 |              0.223 |             0.642 |             6.812 |           6,243.784 |
| Before merge: BioNano Genome Map        |     183 |             0.077 |               29.897 |             33.751 |            61.052 |           175.434 |           6,176.352 |
| Before merge: NGS sequences               |  29,998 |             0.000 |                0.067 |              0.208 |             0.565 |             4.751 |           6,243.792 |
| BNG Genome Map in hybrid scaffold       |     137 |             0.422 |               44.871 |             44.854 |            61.052 |           175.434 |           6,144.991 |
| NGS sequences in hybrid scaffold (CMAP) |  14,321 |             0.028 |                0.272 |              0.405 |             0.619 |             4.751 |           5,797.308 |
| Hybrid scaffold Map                       |     129 |             0.422 |               46.077 |             47.605 |            64.389 |           175.434 |           6,141.024 |
| NGS FASTA sequence in hybrid scaffold   |  14,067 |             0.028 |                0.277 |              0.410 |             0.622 |             4.751 |           5,767.258 |
| Hybrid scaffold FASTA                     |     129 |             0.135 |               46.111 |             47.592 |            65.127 |           177.242 |           6,139.374 |


The `Hybrid scaffold FASTA` was then used for downstream processing.
