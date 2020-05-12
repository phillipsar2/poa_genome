GQ324542.1.fasta
Poa pratensis voucher CAN:Gillespie 6291 internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence
Downloaded from Genbank on 5/11/20

GQ324525.1.fasta
Poa matthewsii voucher OTA:Lloyd s.n. 058920 internal transcribed spacer 1, 5.8S ribosomal RNA gene, and internal transcribed spacer 2, complete sequence
Downloaded from Genbank on 5/12/20

*.consensus.ITS.fasta 
These files are the consensus sequences for the ITS region generated for each induvidual. Reads for each induvidual were mapped to GQ324525.1.fasta, bams were sorted, duplicatees were marked, and read groups were added. SNPs were called with gatk HaplotypeCaller, SNPs were pulled out and hard filtered following GATK best practices. The vcfs were then use to generate a pseduo sequence with gatk FastaAlternateReferenceMaker.

ALL.poa.ITS.fasta
All of the poa ITS consensus sequences in one file.
