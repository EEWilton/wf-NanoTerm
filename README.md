# wf-NanoTerm

NanoTerm is a NextFlow workflow inspired by the tool PhageTerm (Garneau et al, 2017, Scientific Reports).
It is designed to identify the termini of a phage genome based on long-read sequencing results (i.e. Oxford Nanopore).

This workflow uses a set of bioinformatics tools that are all included in the docker image 'wiltone/nanoterm:1.0'.  Alignment of the sequence reads to the reference is done with minimap2.  

Required parameters:
1. --fastq <folder containing processed sequence reads as fastq files>
2. --reference <fasta file of the reference genome>

Optional parameters:
1. --out_dir <diretory to save output files>
2. --seqplat <sequencing platform> (default is nanopore, alternative is illumina)
3. --name <name of the phage to be used in the final report>

This document last updated on June 21, 2023.
Emily Wilton (M.Sc.) - Ph.D. student at University of Manitoba
