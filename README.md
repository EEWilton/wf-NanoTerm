# wf-NanoTerm

NanoTerm is a NextFlow workflow designed to identify the termini of a phage genome based on long-read sequencing results (i.e. Oxford Nanopore).

Inspired by the tool PhageTerm (Garneau et al, 2017, Scientific Reports), NanoTerm calculates a value 'tau', which reflects the percentage of read depth at any position that correspond to the first nucleotide of a sequence read.  A tau value of 1.0 at a specific genome position means that every sequence read covering that position also starts at that position.  A tau value of 0 means that no sequence reads start at that position.  Calculating tau allows for the detection of the physical ends of the sequenced DNA/RNA.

Due to the nature of mapping to a linear reference genome, the first and last positions of the reference genome will have tau values of 1.0.  This is because it is not possible for a read to cover position 1 without also starting at position 1.  To combat this technical artifact, NanoTerm produces five circular permutations of the reference genome, against which the sequence reads are mapped.  By averaging the tau values across the six replicates, these high tau value artifacts are reduced.

![end artifacts](https://github.com/EEWilton/wf-NanoTerm/blob/main/Images/end_artifacts.png) 
Figure 1. An example of tau values mapped across the reference genome, highlighting the artificially high tau values at the first and last positions, compared to the relevantg tau values in the middle indicating the true termini.

![circular permutations](https://github.com/EEWilton/wf-NanoTerm/blob/main/Images/permutations.png)
Figure 2.  An illustration of the five circular permutations of the reference genome compared to the original, with the star indicating the real terminus position that we are looking for.

This workflow uses a set of bioinformatics tools that are all included in the docker image 'wiltone/nanoterm:1.0'.  Alignment of the sequence reads to the reference is done with minimap2.  The samtools package is used to determine total read depth.  The calculations, logic, and final report generation all use in R.  Some data processing was done with Python.

### Usage:

This workflow can be imported to Epi2Me Labs (Oxford Nanopore Techonology) and executed in the Epi2Me Labs GUI.  It can also be run from the command line.

nextflow run /path/to/main.nf --fastq /path/to/fastq/dir --reference /path/to/reference/fasta

Required parameters:
1. --fastq [folder containing processed sequence reads as fastq files]
2. --reference [fasta file of the reference genome]

Optional parameters:
1. --out_dir [directory to save output files]
2. --seqplat [sequencing platform] (default is nanopore, alternative is illumina)
3. --name [name of the phage to be used in the final report] (default is my_phage)

This document last updated on June 21, 2023.

Emily Wilton (M.Sc.) - Ph.D. student at University of Manitoba