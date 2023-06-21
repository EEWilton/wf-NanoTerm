# wf-NanoTerm

NanoTerm is a NextFlow workflow inspired by the tool PhageTerm (Garneau et al, 2017, Scientific Reports).  It is designed to identify the termini of a phage genome based on long-read sequencing results (i.e. Oxford Nanopore).

This workflow uses a set of bioinformatics tools that are all included in the docker image 'wiltone/nanoterm:1.0'.  

Alignment of the sequence reads to the reference is done with minimap2.  The samtools package is used to determine total read depth.  The calculations, logical deductions, and final report generation were all done in R.

Required parameters:
1. --fastq [folder containing processed sequence reads as fastq files]
2. --reference [fasta file of the reference genome]

Optional parameters:
1. --out_dir [diretory to save output files]
2. --seqplat [sequencing platform] (default is nanopore, alternative is illumina)
3. --name [name of the phage to be used in the final report] (default is my_phage)

Usage:
nextflow run /path/to/main.nf --fastq /path/to/fastq/dir --reference /path/to/reference/fasta --name Phage1234

This document last updated on June 21, 2023.

Emily Wilton (M.Sc.) - Ph.D. student at University of Manitoba