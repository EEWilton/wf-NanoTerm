# wf-NanoTerm

NanoTerm is a NextFlow workflow inspired by the tool PhageTerm (Garneau et al, 2017, Scientific Reports).  It is designed to identify the termini of a phage genome based on long-read sequencing results (i.e. Oxford Nanopore).

Specicially, this workflow calculates a value 'tau', which reflects the percentage of read depth at any position that correspond to the first nucleotide of a sequence read.  For example, if a genome position has a read depth of 10 and 5 of those reads start at that position, then the tau value will be 0.5.  A tau value closer to 1 means that most sequence reads that cover the position will start at that position.  Calculating tau allows for the detection of the physical ends of the sequenced DNA/RNA.

This workflow uses a set of bioinformatics tools that are all included in the docker image 'wiltone/nanoterm:1.0'.  Alignment of the sequence reads to the reference is done with minimap2.  The samtools package is used to determine total read depth.  The calculations, logical deductions, and final report generation were all done in R.  Some data processing was done with Python.

# Usage:

This workflow can be imported to Epi2Me Labs (Oxford Nanopore Techonology) and executed in the Epi2Me Labs GUI.  It can also be run from the command line.

nextflow run /path/to/main.nf --fastq /path/to/fastq/dir --reference /path/to/reference/fasta

Required parameters:
1. --fastq [folder containing processed sequence reads as fastq files]
2. --reference [fasta file of the reference genome]

Optional parameters:
1. --out_dir [diretory to save output files]
2. --seqplat [sequencing platform] (default is nanopore, alternative is illumina)
3. --name [name of the phage to be used in the final report] (default is my_phage)

This document last updated on June 21, 2023.

Emily Wilton (M.Sc.) - Ph.D. student at University of Manitoba