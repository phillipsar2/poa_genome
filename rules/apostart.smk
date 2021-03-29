# Identifying the APOSTART region in the reference genome
rule bwa_map:
    input: 
        gene = "data/gene/apostart/APOSTART_sequence.fasta",
        ref = config.ref
    output:
        "data/interm/mapped_bam/apostat.mapped.bam"
    run:
        shell("bwa mem {input.ref} {input.gene} | samtools view -Sb > {output}")
