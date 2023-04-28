#!/usr/bin/env nextflow

ref_ch = Channel.of(params.reference)
seq_ch = Channel.of(params.fastq)
plat_ch = Channel.of(params.seqplat)
name_ch = Channel.of(params.name)
outdir_ch = Channel.of(params.out_dir)

// This process concatenates all of the fastq.gz files in the fastq input directory into one file containing all sequences
process catFastq {
	input:
		path fastq
	
	output:
		path 'all.fastq.gz'
		
	script:
	"""
	cat $fastq/*.fastq.gz > all.fastq.gz
	"""
}

// This process extracts the raw sequence from the fasta file, removing the header
// The output is a text file with just the nucleotide sequence of the reference genome
process rawSeq {
	input:
		path reference
	
	output:
		path 'rawseq.txt'
		
	script:
	"""
	grep -v ">" $reference | tr -d "[:space:]" > rawseq.txt
	"""
}

process refLen {
	input:
		path rawseq
		
	output:
		env genLen
		
	script:
	"""
	genLen=`wc -m < $rawseq`
	"""
}

// This process strips all new lines from the fasta file, except between the header and the sequence
process refSeq {
	input:
		path reference
	
	output:
		path 'refseq.fasta'
		
	script:
	"""
	echo ">input_reference" > refseq.fasta
	grep -v ">" $reference | tr -d "\\n" >> refseq.fasta
	"""
}

// This process produces 5 circular permutations of the reference genome
process permute {
	input:
		path refseq
		
	output:
		path 'circular_permutation1.fasta'
		path 'circular_permutation2.fasta'
		path 'circular_permutation3.fasta'
		path 'circular_permutation4.fasta'
		path 'circular_permutation5.fasta'
	
	script:
	"""
	#!/usr/bin/env python3

	f = open('refseq.fasta', 'r')
	next(f)
	for line in f:
		refseq = str(line.rstrip())
	f.close()

	L = len(refseq)
	br=int(L/6)

	break1 = br
	break2 = br*2
	break3 = br*3
	break4 = br*4
	break5 = br*5

	print(break4)

	seq1 = refseq[break1:] + refseq[:break1]
	seq2 = refseq[break2:] + refseq[:break2]
	seq3 = refseq[break3:] + refseq[:break3]
	seq4 = refseq[break4:] + refseq[:break4]
	seq5 = refseq[break5:] + refseq[:break5]

	sourceFile = open('circular_permutation1.fasta', 'w')
	print('>circular_permutation_1', file = sourceFile)
	print(seq1, file = sourceFile)
	sourceFile.close()

	sourceFile = open('circular_permutation2.fasta', 'w')
	print('>circular_permutation_2', file = sourceFile)
	print(seq2, file = sourceFile)
	sourceFile.close()

	sourceFile = open('circular_permutation3.fasta', 'w')
	print('>circular_permutation_3', file = sourceFile)
	print(seq3, file = sourceFile)
	sourceFile.close()

	sourceFile = open('circular_permutation4.fasta', 'w')
	print('>circular_permutation_4', file = sourceFile)
	print(seq4, file = sourceFile)
	sourceFile.close()

	sourceFile = open('circular_permutation5.fasta', 'w')
	print('>circular_permutation_5', file = sourceFile)
	print(seq5, file = sourceFile)
	sourceFile.close()
	"""
}

// This process maps the input fastq.gz file against the reference genome using minimap2
// This process can be adjusted to suit the sequencing platform used
// The output is a .sam file of the alignment
process mapping {
	input:
		val seqplat
		path refseq
		path all
		path circular_permutation1
		path circular_permutation2
		path circular_permutation3
		path circular_permutation4
		path circular_permutation5

	output:
		path 'aln.sam'
		path 'aln_circ1.sam'
		path 'aln_circ2.sam'
		path 'aln_circ3.sam'
		path 'aln_circ4.sam'
		path 'aln_circ5.sam'
		
	script:
		if( seqplat == 'nanopore' )
			"""
			minimap2 -ax map-ont $refseq $all > aln.sam
			minimap2 -ax map-ont $circular_permutation1 $all > aln_circ1.sam
			minimap2 -ax map-ont $circular_permutation2 $all > aln_circ2.sam
			minimap2 -ax map-ont $circular_permutation3 $all > aln_circ3.sam
			minimap2 -ax map-ont $circular_permutation4 $all > aln_circ4.sam
			minimap2 -ax map-ont $circular_permutation5 $all > aln_circ5.sam
			"""
		
		else if( seqplat == 'illumina' )
			"""
			minimap2 -ax sr $refseq $all > aln.sam
			minimap2 -ax sr $circular_permutation1 $all > aln_circ1.sam
			minimap2 -ax sr $circular_permutation2 $all > aln_circ2.sam
			minimap2 -ax sr $circular_permutation3 $all > aln_circ3.sam
			minimap2 -ax sr $circular_permutation4 $all > aln_circ4.sam
			minimap2 -ax sr $circular_permutation5 $all > aln_circ5.sam
			"""

		else 
			 error "Invalid sequencing platform: $params.seqplat"
}

// This process extracts some statistics from the alignment to the original reference genome
process alignStats {
	input:
		path aln
		path aln_circ1
		path aln_circ2
		path aln_circ3
		path aln_circ4
		path aln_circ5
		
	output:
		env totalReads
		env mappedReads
		env unmappedReads
		env aveReadLen
		env maxReadLen
	
	script:
		"""
		samtools stats $aln | grep ^SN | cut -f 2- > aln_stats.txt 
		
		totalReads=`grep "raw total sequences:" aln_stats.txt | grep -oz [[:digit:]]`
		mappedReads=`grep "reads mapped:" aln_stats.txt | grep -oz [[:digit:]]`
		unmappedReads=`grep "reads unmapped:" aln_stats.txt | grep -oz [[:digit:]]`
		aveReadLen=`grep "average length:" aln_stats.txt | grep -oz [[:digit:]]`
		maxReadLen=`grep "maximum length:" aln_stats.txt | grep -oz [[:digit:]]`
		"""
}

// This process seperates the alignment file by strand, to allow for easier analysis.  It also sorts the .bam files
process strandSep {
	input:
		path aln
		path aln_circ1
		path aln_circ2
		path aln_circ3
		path aln_circ4
		path aln_circ5
		
	output:
		path 'aln_f_sorted.bam'
		path 'aln_r_sorted.bam'
		path 'aln_f_sorted_circ1.bam'
		path 'aln_r_sorted_circ1.bam'
		path 'aln_f_sorted_circ2.bam'
		path 'aln_r_sorted_circ2.bam'
		path 'aln_f_sorted_circ3.bam'
		path 'aln_r_sorted_circ3.bam'
		path 'aln_f_sorted_circ4.bam'
		path 'aln_r_sorted_circ4.bam'
		path 'aln_f_sorted_circ5.bam'
		path 'aln_r_sorted_circ5.bam'

	script:
		"""
		samtools view -F 16 -b $aln > aln_f.bam
		samtools sort aln_f.bam > aln_f_sorted.bam
		
		samtools view -f 16 -b $aln > aln_r.bam
		samtools sort aln_r.bam > aln_r_sorted.bam
		
		samtools view -F 16 -b $aln_circ1 > aln_f_circ1.bam
		samtools sort aln_f_circ1.bam > aln_f_sorted_circ1.bam
		
		samtools view -f 16 -b $aln_circ1 > aln_r_circ1.bam
		samtools sort aln_r_circ1.bam > aln_r_sorted_circ1.bam

		samtools view -F 16 -b $aln_circ2 > aln_f_circ2.bam
		samtools sort aln_f_circ2.bam > aln_f_sorted_circ2.bam
		
		samtools view -f 16 -b $aln_circ2 > aln_r_circ2.bam
		samtools sort aln_r_circ2.bam > aln_r_sorted_circ2.bam
		
		samtools view -F 16 -b $aln_circ3 > aln_f_circ3.bam
		samtools sort aln_f_circ3.bam > aln_f_sorted_circ3.bam
		
		samtools view -f 16 -b $aln_circ3 > aln_r_circ3.bam
		samtools sort aln_r_circ3.bam > aln_r_sorted_circ3.bam
		
		samtools view -F 16 -b $aln_circ4 > aln_f_circ4.bam
		samtools sort aln_f_circ4.bam > aln_f_sorted_circ4.bam
		
		samtools view -f 16 -b $aln_circ4 > aln_r_circ4.bam
		samtools sort aln_r_circ4.bam > aln_r_sorted_circ4.bam
		
		samtools view -F 16 -b $aln_circ5 > aln_f_circ5.bam
		samtools sort aln_f_circ5.bam > aln_f_sorted_circ5.bam
		
		samtools view -f 16 -b $aln_circ5 > aln_r_circ5.bam
		samtools sort aln_r_circ5.bam > aln_r_sorted_circ5.bam
		"""
}

// This process makes a bed file from each sorted bam file, to allow the extraction of starting position coverage
process bed {
	input:
		path aln_f_sorted
		path aln_r_sorted
		path aln_f_sorted_circ1
		path aln_r_sorted_circ1
		path aln_f_sorted_circ2
		path aln_r_sorted_circ2
		path aln_f_sorted_circ3
		path aln_r_sorted_circ3
		path aln_f_sorted_circ4
		path aln_r_sorted_circ4
		path aln_f_sorted_circ5
		path aln_r_sorted_circ5
		
	output:
		path 'aln_f.bed'
		path 'aln_r.bed'
		path 'aln_f_circ1.bed'
		path 'aln_r_circ1.bed'
		path 'aln_f_circ2.bed'
		path 'aln_r_circ2.bed'
		path 'aln_f_circ3.bed'
		path 'aln_r_circ3.bed'
		path 'aln_f_circ4.bed'
		path 'aln_r_circ4.bed'
		path 'aln_f_circ5.bed'
		path 'aln_r_circ5.bed'
		
	script:
	"""
	bamToBed -i $aln_f_sorted > aln_f.bed
	bamToBed -i $aln_r_sorted > aln_r.bed
	
	bamToBed -i $aln_f_sorted_circ1 > aln_f_circ1.bed
	bamToBed -i $aln_r_sorted_circ1 > aln_r_circ1.bed
	
	bamToBed -i $aln_f_sorted_circ2 > aln_f_circ2.bed
	bamToBed -i $aln_r_sorted_circ2 > aln_r_circ2.bed
	
	bamToBed -i $aln_f_sorted_circ3 > aln_f_circ3.bed
	bamToBed -i $aln_r_sorted_circ3 > aln_r_circ3.bed
	
	bamToBed -i $aln_f_sorted_circ4 > aln_f_circ4.bed
	bamToBed -i $aln_r_sorted_circ4 > aln_r_circ4.bed
	
	bamToBed -i $aln_f_sorted_circ5 > aln_f_circ5.bed
	bamToBed -i $aln_r_sorted_circ5 > aln_r_circ5.bed
	"""
}

// This process extracts total sequencing depth across the genome
process cov {	
	input:
		path aln_f_sorted
		path aln_r_sorted
		path aln_f_sorted_circ1
		path aln_r_sorted_circ1
		path aln_f_sorted_circ2
		path aln_r_sorted_circ2
		path aln_f_sorted_circ3
		path aln_r_sorted_circ3
		path aln_f_sorted_circ4
		path aln_r_sorted_circ4
		path aln_f_sorted_circ5
		path aln_r_sorted_circ5
		
	output:
		path 'aln_f_cov.txt'
		path 'aln_r_cov.txt'
		path 'aln_f_cov_circ1.txt'
		path 'aln_r_cov_circ1.txt'
		path 'aln_f_cov_circ2.txt'
		path 'aln_r_cov_circ2.txt'
		path 'aln_f_cov_circ3.txt'
		path 'aln_r_cov_circ3.txt'
		path 'aln_f_cov_circ4.txt'
		path 'aln_r_cov_circ4.txt'
		path 'aln_f_cov_circ5.txt'
		path 'aln_r_cov_circ5.txt'
	
	script:
	"""
	samtools depth -a -d 0 -o aln_f_cov.txt $aln_f_sorted
	samtools depth -a -d 0 -o aln_r_cov.txt $aln_r_sorted
	
	samtools depth -a -d 0 -o aln_f_cov_circ1.txt $aln_f_sorted_circ1
	samtools depth -a -d 0 -o aln_r_cov_circ1.txt $aln_r_sorted_circ1
	
	samtools depth -a -d 0 -o aln_f_cov_circ2.txt $aln_f_sorted_circ2
	samtools depth -a -d 0 -o aln_r_cov_circ2.txt $aln_r_sorted_circ2
	
	samtools depth -a -d 0 -o aln_f_cov_circ3.txt $aln_f_sorted_circ3
	samtools depth -a -d 0 -o aln_r_cov_circ3.txt $aln_r_sorted_circ3
	
	samtools depth -a -d 0 -o aln_f_cov_circ4.txt $aln_f_sorted_circ4
	samtools depth -a -d 0 -o aln_r_cov_circ4.txt $aln_r_sorted_circ4
	
	samtools depth -a -d 0 -o aln_f_cov_circ5.txt $aln_f_sorted_circ5
	samtools depth -a -d 0 -o aln_r_cov_circ5.txt $aln_r_sorted_circ5
	"""
}

// This process calculates tau for each position in each circular permutation
process tau {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true
	
	input:
		val genLen
		path aln_f
		path aln_r
		path aln_f_circ1
		path aln_r_circ1
		path aln_f_circ2
		path aln_r_circ2
		path aln_f_circ3
		path aln_r_circ3
		path aln_f_circ4
		path aln_r_circ4
		path aln_f_circ5
		path aln_r_circ5
		path aln_f_cov
		path aln_r_cov
		path aln_f_cov_circ1
		path aln_r_cov_circ1
		path aln_f_cov_circ2
		path aln_r_cov_circ2
		path aln_f_cov_circ3
		path aln_r_cov_circ3
		path aln_f_cov_circ4
		path aln_r_cov_circ4
		path aln_f_cov_circ5
		path aln_r_cov_circ5
		
	output:
		path 'tau.csv'
		path 'tau_circ1.csv'
		path 'tau_circ2.csv'
		path 'tau_circ3.csv'
		path 'tau_circ4.csv'
		path 'tau_circ5.csv'
		path 'all_tau.csv'
		path 'tau_stats.csv'
		
		
	script:
	"""
	#!/usr/bin/env Rscript
	
	library("tidyverse")
	library("dplyr")

	len <- $genLen

	break1 <- as.integer(len / 6)
	break2 <- break1 * 2
	break3 <- break1 * 3
	break4 <- break1 * 4
	break5 <- break1 * 5
		
	pos <- seq(from = 1,to = len, by = 1)

	f_bed <- read.table("$aln_f", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	f_bed <- setNames(f_bed, c("chr","start","end"))
	f_bed['start'] <- f_bed['start'] + 1

	r_bed <- read.table("$aln_r", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	r_bed <- setNames(r_bed, c("chr","start","end"))
	r_bed['start'] <- r_bed['start'] + 1

	f_bed_circ1 <- read.table("$aln_f_circ1", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	f_bed_circ1 <- setNames(f_bed_circ1, c("chr","start","end"))
	f_bed_circ1['start'] <- f_bed_circ1['start'] + 1

	r_bed_circ1 <- read.table("$aln_r_circ1", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	r_bed_circ1 <- setNames(r_bed_circ1, c("chr","start","end"))
	r_bed_circ1['start'] <- r_bed_circ1['start'] + 1

	f_bed_circ2 <- read.table("$aln_f_circ2", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	f_bed_circ2 <- setNames(f_bed_circ2, c("chr","start","end"))
	f_bed_circ2['start'] <- f_bed_circ2['start'] + 1

	r_bed_circ2 <- read.table("$aln_r_circ2", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	r_bed_circ2 <- setNames(r_bed_circ2, c("chr","start","end"))
	r_bed_circ2['start'] <- r_bed_circ2['start'] + 1

	f_bed_circ3 <- read.table("$aln_f_circ3", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	f_bed_circ3 <- setNames(f_bed_circ3, c("chr","start","end"))
	f_bed_circ3['start'] <- f_bed_circ3['start'] + 1

	r_bed_circ3 <- read.table("$aln_r_circ3", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	r_bed_circ3 <- setNames(r_bed_circ3, c("chr","start","end"))
	r_bed_circ3['start'] <- r_bed_circ3['start'] + 1
	
	f_bed_circ4 <- read.table("$aln_f_circ4", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	f_bed_circ4 <- setNames(f_bed_circ4, c("chr","start","end"))
	f_bed_circ4['start'] <- f_bed_circ4['start'] + 1

	r_bed_circ4 <- read.table("$aln_r_circ4", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	r_bed_circ4 <- setNames(r_bed_circ4, c("chr","start","end"))
	r_bed_circ4['start'] <- r_bed_circ4['start'] + 1
	
	f_bed_circ5 <- read.table("$aln_f_circ5", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	f_bed_circ5 <- setNames(f_bed_circ5, c("chr","start","end"))
	f_bed_circ5['start'] <- f_bed_circ5['start'] + 1

	r_bed_circ5 <- read.table("$aln_r_circ5", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	r_bed_circ5 <- setNames(r_bed_circ5, c("chr","start","end"))
	r_bed_circ5['start'] <- r_bed_circ5['start'] + 1


	f_cov <- read.table("$aln_f_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	f_cov <- setNames(f_cov, c("chr","pos","cov"))
	f_cov['strand'] <- "f"

	r_cov <- read.table("$aln_r_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	r_cov <- setNames(r_cov, c("chr","pos","cov"))
	r_cov['strand'] <- "r"
	
	f_cov_circ1 <- read.table("$aln_f_cov_circ1", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	f_cov_circ1 <- setNames(f_cov_circ1, c("chr","pos","cov"))
	f_cov_circ1['strand'] <- "f"

	r_cov_circ1 <- read.table("$aln_r_cov_circ1", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	r_cov_circ1 <- setNames(r_cov_circ1, c("chr","pos","cov"))
	r_cov_circ1['strand'] <- "r"
	
	f_cov_circ2 <- read.table("$aln_f_cov_circ2", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	f_cov_circ2 <- setNames(f_cov_circ2, c("chr","pos","cov"))
	f_cov_circ2['strand'] <- "f"

	r_cov_circ2 <- read.table("$aln_r_cov_circ2", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	r_cov_circ2 <- setNames(r_cov_circ2, c("chr","pos","cov"))
	r_cov_circ2['strand'] <- "r"

	f_cov_circ3 <- read.table("$aln_f_cov_circ3", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	f_cov_circ3 <- setNames(f_cov_circ3, c("chr","pos","cov"))
	f_cov_circ3['strand'] <- "f"

	r_cov_circ3 <- read.table("$aln_r_cov_circ3", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	r_cov_circ3 <- setNames(r_cov_circ3, c("chr","pos","cov"))
	r_cov_circ3['strand'] <- "r"

	f_cov_circ4 <- read.table("$aln_f_cov_circ4", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	f_cov_circ4 <- setNames(f_cov_circ4, c("chr","pos","cov"))
	f_cov_circ4['strand'] <- "f"

	r_cov_circ4 <- read.table("$aln_r_cov_circ4", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	r_cov_circ4 <- setNames(r_cov_circ4, c("chr","pos","cov"))
	r_cov_circ4['strand'] <- "r"	

	f_cov_circ5 <- read.table("$aln_f_cov_circ5", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	f_cov_circ5 <- setNames(f_cov_circ5, c("chr","pos","cov"))
	f_cov_circ5['strand'] <- "f"

	r_cov_circ5 <- read.table("$aln_r_cov_circ5", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
	r_cov_circ5 <- setNames(r_cov_circ5, c("chr","pos","cov"))
	r_cov_circ5['strand'] <- "r"	
		
	
	f_spc <- data.frame(pos = pos, SPC = 0)
	r_spc <- data.frame(pos = pos, SPC = 0)

	for (nt in pos){
		f_spc[nt,'SPC'] <- length(which(f_bed['start']==nt))
		}
	
	for (nt in pos){
		r_spc[nt,'SPC'] <- length(which(r_bed['end']==nt))
		}
	
	f_tau <- merge(f_cov, f_spc, by="pos")
	r_tau <- merge(r_cov, r_spc, by="pos")
	
	tau <- rbind(f_tau, r_tau)
	tau['tau'] <- tau['SPC'] / tau['cov']
	
	for (i in 1:dim(tau)[1]) {
		tau[i, 'pos_fix'] <- tau[i, 'pos'] 
	}

	for (i in 1:dim(tau)[1]) {
		if (tau[i, 'pos_fix'] > len) {
			tau[i, 'pos_adj'] <- tau[i, 'pos_fix'] - len
		} else {
			tau[i, 'pos_adj'] <- tau[i, 'pos_fix']
		}
	}
	
	tau <- subset(tau, select=-pos_fix)
	write.csv(tau, "tau.csv")
	
	f_spc_circ1 <- data.frame(pos = pos, SPC = 0)
	r_spc_circ1 <- data.frame(pos = pos, SPC = 0)

	for (nt in pos){
		f_spc_circ1[nt,'SPC'] <- length(which(f_bed_circ1['start']==nt))
		}
	
	for (nt in pos){
		r_spc_circ1[nt,'SPC'] <- length(which(r_bed_circ1['end']==nt))
		}
	
	f_tau_circ1 <- merge(f_cov_circ1, f_spc_circ1, by="pos")
	r_tau_circ1 <- merge(r_cov_circ1, r_spc_circ1, by="pos")
	
	tau_circ1 <- rbind(f_tau_circ1, r_tau_circ1)
	tau_circ1['tau'] <- tau_circ1['SPC'] / tau_circ1['cov']
	
	for (i in 1:dim(tau_circ1)[1]) {
		tau_circ1[i, 'pos_fix'] <- tau_circ1[i, 'pos'] + break1 
	}

	for (i in 1:dim(tau_circ1)[1]) {
		if (tau_circ1[i, 'pos_fix'] > len) {
				tau_circ1[i, 'pos_adj'] <- tau_circ1[i, 'pos_fix'] - len
		} else {
			tau_circ1[i, 'pos_adj'] <- tau_circ1[i, 'pos_fix']
		}
	}
	
	tau_circ1 <- subset(tau_circ1, select=-pos_fix)
	write.csv(tau_circ1, "tau_circ1.csv")
	
	f_spc_circ2 <- data.frame(pos = pos, SPC = 0)
	r_spc_circ2 <- data.frame(pos = pos, SPC = 0)

	for (nt in pos){
		f_spc_circ2[nt,'SPC'] <- length(which(f_bed_circ2['start']==nt))
		}
	
	for (nt in pos){
		r_spc_circ2[nt,'SPC'] <- length(which(r_bed_circ2['end']==nt))
		}
	
	f_tau_circ2 <- merge(f_cov_circ2, f_spc_circ2, by="pos")
	r_tau_circ2 <- merge(r_cov_circ2, r_spc_circ2, by="pos")
	
	tau_circ2 <- rbind(f_tau_circ2, r_tau_circ2)
	tau_circ2['tau'] <- tau_circ2['SPC'] / tau_circ2['cov']

	for (i in 1:dim(tau_circ2)[1]) {
		tau_circ2[i, 'pos_fix']<- tau_circ2[i, 'pos'] + break2 
	}

	for (i in 1:dim(tau_circ2)[1]) {
		if (tau_circ2[i, 'pos_fix'] > len) {
			tau_circ2[i, 'pos_adj'] <- tau_circ2[i, 'pos_fix'] - len
		} else {
		tau_circ2[i, 'pos_adj'] <- tau_circ2[i, 'pos_fix']
		}
	}

	tau_circ2 <- subset(tau_circ2, select=-pos_fix)
	write.csv(tau_circ2, "tau_circ2.csv")

	f_spc_circ3 <- data.frame(pos = pos, SPC = 0)
	r_spc_circ3 <- data.frame(pos = pos, SPC = 0)

	for (nt in pos){
		f_spc_circ3[nt,'SPC'] <- length(which(f_bed_circ3['start']==nt))
		}
	
	for (nt in pos){
		r_spc_circ3[nt,'SPC'] <- length(which(r_bed_circ3['end']==nt))
		}
	
	f_tau_circ3 <- merge(f_cov_circ3, f_spc_circ3, by="pos")
	r_tau_circ3 <- merge(r_cov_circ3, r_spc_circ3, by="pos")
	
	tau_circ3 <- rbind(f_tau_circ3, r_tau_circ3)
	tau_circ3['tau'] <- tau_circ3['SPC'] / tau_circ3['cov']
	
	for (i in 1:dim(tau_circ3)[1]) {
		tau_circ3[i, 'pos_fix'] <- tau_circ3[i, 'pos'] + break3 
	}

	for (i in 1:dim(tau_circ3)[1]) {
		if (tau_circ3[i, 'pos_fix'] > len) {
			tau_circ3[i, 'pos_adj'] <- tau_circ3[i, 'pos_fix'] - len
		} else {
			tau_circ3[i, 'pos_adj'] <- tau_circ3[i, 'pos_fix']
		}
	}
	
	tau_circ3 <- subset(tau_circ3, select=-pos_fix)
	write.csv(tau_circ3, "tau_circ3.csv")

	f_spc_circ4 <- data.frame(pos = pos, SPC = 0)
	r_spc_circ4 <- data.frame(pos = pos, SPC = 0)

	for (nt in pos){
		f_spc_circ4[nt,'SPC'] <- length(which(f_bed_circ4['start']==nt))
		}
	
	for (nt in pos){
		r_spc_circ4[nt,'SPC'] <- length(which(r_bed_circ4['end']==nt))
		}
	
	f_tau_circ4 <- merge(f_cov_circ4, f_spc_circ4, by="pos")
	r_tau_circ4 <- merge(r_cov_circ4, r_spc_circ4, by="pos")
	
	tau_circ4 <- rbind(f_tau_circ4, r_tau_circ4)
	tau_circ4['tau'] <- tau_circ4['SPC'] / tau_circ4['cov']
	
	for (i in 1:dim(tau_circ4)[1]) {
		tau_circ4[i, 'pos_fix'] <- tau_circ4[i, 'pos'] + break4 
	}

	for (i in 1:dim(tau_circ4)[1]) {
		if (tau_circ4[i, 'pos_fix'] > len) {
			tau_circ4[i, 'pos_adj'] <- tau_circ4[i, 'pos_fix'] - len
		} else {
			tau_circ4[i, 'pos_adj'] <- tau_circ4[i, 'pos_fix']
		}
	}
	
	tau_circ4 <- subset(tau_circ4, select=-pos_fix)
	write.csv(tau_circ4, "tau_circ4.csv")
	
	f_spc_circ5 <- data.frame(pos = pos, SPC = 0)
	r_spc_circ5 <- data.frame(pos = pos, SPC = 0)

	for (nt in pos){
		f_spc_circ5[nt,'SPC'] <- length(which(f_bed_circ5['start']==nt))
		}
	
	for (nt in pos){
		r_spc_circ5[nt,'SPC'] <- length(which(r_bed_circ5['end']==nt))
		}
	
	f_tau_circ5 <- merge(f_cov_circ5, f_spc_circ5, by="pos")
	r_tau_circ5 <- merge(r_cov_circ5, r_spc_circ5, by="pos")
	
	tau_circ5 <- rbind(f_tau_circ5, r_tau_circ5)
	tau_circ5['tau'] <- tau_circ5['SPC'] / tau_circ5['cov']
	
	for (i in 1:dim(tau_circ5)[1]) {
		tau_circ5[i, 'pos_fix'] <- tau_circ5[i, 'pos'] + break5 
	}

	for (i in 1:dim(tau_circ5)[1]) {
		if (tau_circ5[i, 'pos_fix'] > len) {
			tau_circ5[i, 'pos_adj'] <- tau_circ5[i, 'pos_fix'] - len
		} else {
			tau_circ5[i, 'pos_adj'] <- tau_circ5[i, 'pos_fix']
		}
	}
	
	tau_circ5 <- subset(tau_circ5, select=-pos_fix)
	write.csv(tau_circ5, "tau_circ5.csv")
	
	all_tau <- rbind(tau, tau_circ1, tau_circ2, tau_circ3, tau_circ4, tau_circ5)
		
	for (i in 1:dim(all_tau)[1]) {
		if (all_tau[i, 'strand'] == "f" & all_tau[i, 'pos'] == 1 & all_tau[i, 'tau'] == 1){
			all_tau[i, 'tau'] <- NA
		} else if (all_tau[i, 'strand']  == "r" & all_tau[i, 'pos'] == len & all_tau[i, 'tau'] == 1){
			all_tau[i, 'tau'] <- NA
		}
	}
	
	write.csv(all_tau, "all_tau.csv")
	
	tau_stats <- all_tau %>% 
		group_by(pos_adj, strand) %>%
		summarise(
			avg_tau = mean(tau, na.rm=TRUE),
			sd = sd(tau, na.rm=TRUE))
	write.csv(tau_stats, "tau_stats.csv")
	"""	
}

process classify {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true

	input:
		val genLen
		path tau
		path tau_circ1
		path tau_circ2
		path tau_circ3
		path tau_circ4
		path tau_circ5
		path all_tau
		path tau_stats
		val name
		val totalReads
		val mappedReads
		val unmappedReads
		val aveReadLen
		val maxReadLen
		path fastq
		path reference
		val seqplat

	output:
		path 'classification.csv'
		
	script:
	"""
	#!/usr/bin/env Rscript

	library("tidyverse")
	library("dplyr")

	len <- $genLen

	tau_stats <- read.csv("$tau_stats")
	tau_stats <- subset(tau_stats, select=-X)

	top_sig_plus <- tau_stats %>%
  	filter(strand == "f" & avg_tau >= 0.1) %>%
  		arrange(desc(avg_tau))

	top_sig_minus <- tau_stats %>%
  		filter(strand == "r" & avg_tau >= 0.1) %>%
  		arrange(desc(avg_tau))

	num_sig_plus <- nrow(top_sig_plus)
	num_sig_minus <- nrow(top_sig_minus)

	if (num_sig_plus == 0){
  		plus_term = NA
  		plus_term_tau = NA
	} else {
  		plus_term <- as.integer(top_sig_plus[1,1])
  		plus_term_tau <- top_sig_plus[1,3]
	}

	if (num_sig_minus == 0){
 		 minus_term = NA
 	 	minus_term_tau = NA
	} else {
 		minus_term <- as.integer(top_sig_minus[1,1])
  		minus_term_tau <- top_sig_minus[1,3]
	}

	if (num_sig_plus > 1 | num_sig_minus > 1) {
 		 peaks = "multiple"
	} else if (num_sig_plus == 1 & num_sig_minus == 1) {
  		peaks = "two"
	} else if ((num_sig_plus == 0 | num_sig_minus == 0) & (num_sig_plus == 1 | num_sig_minus == 1)) {
  		peaks = "one"
	} else if (num_sig_plus == 0 & num_sig_minus == 0) {
 		 peaks = "none"
	} else {
  		peaks = NA
	}

	if (peaks == "multiple"){
  		if (plus_term_tau < 0.35 & minus_term_tau < 0.35) {
    		peaks = "multiple"
  		} else {
    		peaks = "two"
  		}
	}

	if (peaks == "two") {
  		if (plus_term == 1 | minus_term == len){
  			preclass = "terminal"
  		} else if (plus_term > 1 & minus_term < len) {
  			preclass = "internal"
  		}
	}

	if (peaks == "two" & preclass == "internal") {
  		term_dist <- abs(minus_term - plus_term)
  		print(term_dist)
  		if (term_dist <= 20) {
    		class = "COS"
    		print(class)
    	if (plus_term < minus_term) {
      		subclass = "5 prime"
   		 } else {
     		subclass = "3 prime"
    	}
  		 print(subclass) 
  	} else {
    	class = "DTR"
   		print(class)
   		if (term_dist < 1000) {
      		subclass = "short"
    	} else {
     		 subclass = "long"
    	}
    	print(subclass)  
    	}
	}

	if (peaks == "one") {
  		if (is.na(plus_term)){
    		if (minus_term == len | minus_term == 1) {
     			preclass = "terminal"
    		} else {
      			preclass = "internal"
    		}
  		} else if (is.na(minus_term)) {
    		if (plus_term == 1 | plus_term == len) {
      			preclass = "terminal"
    		} else {
      		preclass = "internal"
    		}
 		}
	}


	if (peaks == "none"){

	}

	if (peaks == "multiple") {

	}

	category <- c("len", "num_sig_plus", "num_sig_minus", "peaks",
                  "plus_term", "plus_term_tau",
                  "minus_term", "minus_term_tau",
                  "preclass", "subclass", "term_dist")
    value <- c(len, num_sig_plus, num_sig_minus, peaks, 
                plus_term, plus_term_tau,
                minus_term, minus_term_tau,
                preclass, subclass, term_dist)

	# Convert nested list to the dataframe by columns
	df <- data.frame(category, value)
	df

	write.csv(df, "classification.csv")

	"""
}

process terminalReads {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true
	
	input:
		path aln
		path aln_circ1
		path aln_circ2
		path aln_circ3
		path aln_circ4
		path aln_circ5
		path classification
		
	output:
		path 'term_aln.bam'

	shell:
	'''
	samtools view -b !{aln} > aln.bam
	
	samtools sort aln.bam > aln_sorted.bam
	samtools index aln_sorted.bam
	
	plus_term="$(cat !{classification} | tr -d '"' | awk -F "," '$2 == "plus_term" {print $3}')"
	minus_term="$(cat !{classification} | tr -d '"' | awk -F "," '$2 == "minus_term" {print $3}')"

	echo $plus_term
	echo $minus_term

	region="$(echo "input_reference:${plus_term}-${minus_term}")"

	echo $region

	samtools view -b aln_sorted.bam $region > term_aln.bam
	'''
}

process report {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true

	input:
		
	output:
		path 'report.docx'

	script:
	"""
	#!/usr/bin/env Rscript

	library("ggplot2")
	library("tidyverse")
	library("dplyr")
	library("ggrepel")
	library("ggthemes")
	library("scales")
	library("officer")
	library("mschart")
	library("zoo")
	library("stringr")
	library("seqinr")
	"""
}
	

process doc2pdf {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true
	
	input:
		path report
		
	output:
		path 'report.pdf', optional: true
	
	script:
	"""
	libreoffice --headless --convert-to pdf $report --outdir $params.out_dir
	"""	
}

workflow {
	allseq_ch = catFastq(seq_ch)
	raw_ch = rawSeq(ref_ch)
	len_ch = refLen(raw_ch)
	refseq_ch = refSeq(ref_ch)
	circ_ch = permute(refseq_ch)
	aln_ch = mapping(plat_ch, ref_ch, allseq_ch, circ_ch)
	alnstats_ch = alignStats(aln_ch)
	sep_ch = strandSep(aln_ch)
	bed_ch = bed(sep_ch)
	cov_ch = cov(sep_ch)
	tau_ch = tau(len_ch, bed_ch, cov_ch)
	class_ch = classify(len_ch, tau_ch, name_ch, alnstats_ch, seq_ch, ref_ch, plat_ch)
	termreads_ch = terminalReads(aln_ch, class_ch)
	
}