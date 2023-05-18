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

// This process returns the length of the reference genome
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
// It also outputs the file name of the reference file for use in the report
process alignStats {
	input:
		path aln
		path aln_circ1
		path aln_circ2
		path aln_circ3
		path aln_circ4
		path aln_circ5
		path reference
		
	output:
		env totalReads
		env mappedReads
		env unmappedReads
		env aveReadLen
		env maxReadLen
		env refName
	
	script:
		"""
		samtools stats $aln | grep ^SN | cut -f 2- > aln_stats.txt 
		
		totalReads=`grep "raw total sequences:" aln_stats.txt | grep -oz [[:digit:]]`
		mappedReads=`grep "reads mapped:" aln_stats.txt | grep -oz [[:digit:]]`
		unmappedReads=`grep "reads unmapped:" aln_stats.txt | grep -oz [[:digit:]]`
		aveReadLen=`grep "average length:" aln_stats.txt | grep -oz [[:digit:]]`
		maxReadLen=`grep "maximum length:" aln_stats.txt | grep -oz [[:digit:]]`

		refName=`basename $reference`
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
	
	write.csv(all_tau, "all_tau.csv")
	
	for (i in 1:dim(all_tau)[1]) {
		if (!is.na(all_tau[i, 'tau'])) {
			if (all_tau[i, 'strand'] == "f" & all_tau[i, 'pos'] == 1 & all_tau[i, 'tau'] == 1){
				all_tau[i, 'tau'] <- NA
			} else if (all_tau[i, 'strand']  == "r" & all_tau[i, 'pos'] == len & all_tau[i, 'tau'] == 1){
				all_tau[i, 'tau'] <- NA
			}
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

// This process classifies the phage based on the termini, distance, etc
process classify {
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
		val refName
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

	print(len)

	tau_stats <- read.csv("$tau_stats")
	tau_stats <- subset(tau_stats, select=-X)

	tau_circ3 <- read.csv("$tau_circ3")
	tau_circ3 <- subset(tau_circ3, select=-X)

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
  		if (plus_term <= 5 | minus_term >= (len - 5) | minus_term <= 5 | plus_term >= (len - 5)){
  			location = "terminal"
			subclass = "NA"
			term_dist = "NA"
			class = "NA"
  		} else if (plus_term > 5 & minus_term < (len - 5)) {
  			location = "internal"
  		}
	}

		if (peaks == "one") {
  		if (is.na(plus_term)){
    		if (plus_term <= 5 | minus_term >= (len - 5) | minus_term <= 5 | plus_term >= (len - 5)) {
     			location = "terminal"
				subclass = "NA"
    		} else {
      			location = "internal"
    		}
  		} else if (is.na(minus_term)) {
    		if (plus_term <= 5 | minus_term >= (len - 5) | minus_term <= 5 | plus_term >= (len - 5)) {
      			location = "terminal"
				subclass = "NA"
    		} else {
      		location = "internal"
    		}
 		}
	}

	if (peaks == "none"){
		location = NA
		class = "headful without packaging site or Mu-like"
		subclass = NA
		term_dist = NA
	}

	if (peaks == "multiple") {
		location = NA
		class = "multiple"
		subclass = NA
		term_dist = NA
	}

	if (peaks == "two" & location == "internal") {
		term_dist <- abs(minus_term - plus_term)
  		if (term_dist <= 20) {
    		class = "COS"
    		if (plus_term < minus_term) {
      			subclass = "5 prime"
   			} else {
      			subclass = "3 prime"
   			}
  		} else {
   			class = "DTR"
   		    if (term_dist < 1000) {
      			subclass = "short"
    		} else {
      			subclass = "long"
    		}
	    }
	}

	if (peaks == "two" & location == "terminal") {
  		term_dist_ext <- len - (abs(minus_term - plus_term))
		term_dist_int <- abs(minus_term - plus_term)
		dists <- c(term_dist_ext, term_dist_int)
		term_dist <- min(dists)
 		if (term_dist <= 20) {
    		class = "COS"
    	if (plus_term < minus_term) {
      		subclass = "5 prime"
   		 } else {
     		subclass = "3 prime"
    	}
  	} else {
    	class = "DTR"
   		if (term_dist < 1000) {
      		subclass = "short"
    	} else {
     		 subclass = "long"
    	} 
    	}
	}

	if (peaks == "one" & location == "internal") {
		term_dist = NA
		class = "pac"
  		if (is.na(plus_term)) {
			subclass = "reverse"
		} else {
			subclass = "forward"
		}
	}

	if (peaks == "one" & location == "terminal") {
  		term_dist = NA
 		class = "pac"
  		if (is.na(plus_term)) {
			subclass = "reverse"
		} else {
			subclass = "forward"
		}
	}

	plus_term_circ3 <- tau_circ3 %>%
						filter(pos_adj == plus_term & strand == "f") %>%
						select(pos)
	plus_term_circ3 <- as.integer(plus_term_circ3)

	minus_term_circ3 <- tau_circ3 %>%
						filter(pos_adj == minus_term & strand == "r") %>%
						select(pos)
	minus_term_circ3 <- as.integer(minus_term_circ3)

	classification <- data.frame(len, num_sig_plus, num_sig_minus, peaks, 
                             plus_term, plus_term_tau,
                             minus_term, minus_term_tau,
                             location, subclass, class, term_dist,
							 plus_term_circ3, minus_term_circ3)

	write.csv(classification, "classification.csv")
	"""
}

// This process outputs the reads that align with the terminal region
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
		path 'term_aln_circ3.bam'
		env terminalReads
		
	shell:
	'''
	location="$(cat !{classification} | tr -d '"' | awk -F "," '$1 == "1" {print $10}')"
	peaks="$(cat !{classification} | tr -d '"' | awk -F "," '$1 == "1" {print $5}')"
	plus_term="$(cat !{classification} | tr -d '"' | awk -F "," '$1 == "1" {print $6}')"
	minus_term="$(cat !{classification} | tr -d '"' | awk -F "," '$1 == "1" {print $8}')"
	plus_term_circ3="$(cat !{classification} | tr -d '"' | awk -F "," '$1 == "1" {print $14}')"
	minus_term_circ3="$(cat !{classification} | tr -d '"' | awk -F "," '$1 == "1" {print $15}')"

	region=NA
	terminalReads=FALSE
	if [ $location == "internal" ]
	then
		if [ $peaks == "two" ]
		then
			terminalReads=TRUE
			if [ $plus_term -lt $minus_term ]
			then
				region="$(echo input_reference:$plus_term-$minus_term)"
			elif [ $plus_term -gt $minus_term ]
			then
				region="$(echo input_reference:$minus_term-$plus_term)"
			fi
		fi
		if [ $peaks == "one" ]
		then
			terminalReads=TRUE
			if [ plus_term == NA ]
			then
				region="$(echo input_reference:($minus_term - 10)-($minus_term + 10))"
			elif [ minus_term == NA ]
			then
				region="$(echo input_reference:($plus_term - 10)-($plus_term + 10))"
			fi
		fi
		samtools view -b !{aln} > aln.bam
		samtools sort aln.bam > aln_sorted.bam
		samtools index aln_sorted.bam
		samtools view -b aln_sorted.bam \"$region\" > term_aln.bam
		touch term_aln_circ3.bam
	fi

	if [ $location == "terminal" ]
	then
		if [ $peaks == "two" ]
		then
			terminalReads=TRUE
			if [ $plus_term_circ3 -lt $minus_term_circ3 ]
			then
				region="$(echo circular_permutation_3:$plus_term_circ3-$minus_term_circ3)"
			elif [ $plus_term_circ3 -gt $minus_term_circ3 ]
			then
				region="$(echo circular_permutation_3:$minus_term_circ3-$plus_term_circ3)"
			fi
		fi
		if [ $peaks == "one" ]
		then
			terminalReads=TRUE
			if [ $plus_term_circ3 == NA ]
			then
				region="$(echo circular_permutation_3:($minus_term_circ3 - 10)-($minus_term_circ3 + 10))"
			elif [ $minus_term_circ3 == NA ]
			then
				region="$(echo circular_permutation_3:($plus_term_circ3 - 10)-($plus_term_circ3 + 10))"
			fi
		fi
		samtools view -b !{aln_circ3} > aln_circ3.bam
		samtools sort aln_circ3.bam > aln_circ3_sorted.bam
		samtools index aln_circ3_sorted.bam
		samtools view -b aln_circ3_sorted.bam \"$region\" > term_aln_circ3.bam
		touch term_aln.bam
	fi

	if [ $location == "internal" ]
	then
		touch term_aln.bam
		touch term_aln_circ3.bam
	fi
	'''
}

// This process outputs a .docx report on the findings
process report {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true

	input:
		path reference
		val name
		val seqplat
		val totalReads
		val mappedReads
		val unmappedReads
		val aveReadLen
		val maxReadLen
		val refName
		path classification
		path tau
		path tau_circ1
		path tau_circ2
		path tau_circ3
		path tau_circ4
		path tau_circ5
		path all_tau
		path tau_stats
		path term_aln
		path term_aln_circ3
		val terminalReads

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
	library("GenomicRanges")
	library("GenomicAlignments")
	library("ggbio")

	date <- format(Sys.Date())
	name <- as.character("$name")
	totalReads <- as.integer("$totalReads")
	mappedReads <- as.integer("$mappedReads")
	unmappedReads <- as.integer("$unmappedReads")
	aveReadLen <- as.integer("$aveReadLen")
	maxReadLen <- as.integer("$maxReadLen")
	percentMapped <- format(((mappedReads / totalReads) * 100), digits=3)
	percentUnmapped <- format(((unmappedReads / totalReads) * 100), digits=3)

	tau <- read.csv("$tau")
	tau <- subset(tau, select=-X)

	tau_stats <- read.csv("$tau_stats")
	tau_stats <- subset(tau_stats, select=-X)

	classification <- read.csv("$classification")
	classification <- subset(classification, select=-X)
	
	terminalReads <- $terminalReads

	len <- classification[1,1]
	num_sig_plus <- classification[1,2]
	num_sig_minus <- classification[1,3]
	peaks <- classification[1,4]
	plus_term <- classification[1,5]
	plus_term_tau <- classification[1,6]
	minus_term <- classification[1,7]
	minus_term_tau <- classification[1,8]
	location <- classification[1,9]
	subclass <- classification[1,10]
	class <- classification[1,11]
	term_dist <- classification[1,12]
	plus_term_circ3 <- classification[1,13]
	minus_term_circ3 <- classification[1,14]
	
	top_5_plus <- tau_stats %>% 
  		filter(strand == "f") %>%
 		arrange(desc(avg_tau)) %>%
  		head(5)

	top_5_minus <- tau_stats %>% 
  		filter(strand == "r") %>%
  		arrange(desc(avg_tau)) %>%
  		head(5)

	table <- rbind(top_5_plus, top_5_minus)
	colnames(table) <- c('Position','Strand','Average Tau', 'Standard Deviation')
	table[table=="f"] <- "plus"
	table[table=="r"] <- "minus"
	table <- format(table, digits = 2)

	colours <- c("plus" = "#1B9E77", "minus" = "#7570B3", "both" = "brown4")
	window <- 0.01 * len

	sum_cov <- tau %>%
            group_by(pos) %>%
            summarise(covs = sum(cov))

	if (terminalReads == TRUE & location == "internal") {
		strand.labs <- c("Strand: plus", "Strand: minus")
		names(strand.labs) <- c("+", "-")
		term_aln <- readGAlignments("$term_aln", use.names = TRUE)
		ggterm <- autoplot(term_aln, xlab="Reference genome position",
			geom = "rect", aes(fill = strand, colour = strand), show.legend = FALSE) +
 			facet_grid(strand ~ ., labeller = labeller(strand = strand.labs)) +
  			theme_calc() + 
			geom_vline(xintercept=plus_term, linetype="dashed", colour="springgreen3", linewidth=1) + 
  			geom_vline(xintercept=minus_term, linetype="dashed", colour="violet", linewidth=1)
		png("ggterm.png", width = 6, height = 8, units = "in", res = 300)
		print(ggterm)
		dev.off()
	}

	if (terminalReads == TRUE & location == "terminal") {
		term_aln <- readGAlignments("$term_aln_circ3", use.names = TRUE)
		strand.labs <- c("Strand: plus", "Strand: minus")
		names(strand.labs) <- c("+", "-")
		ggterm <-autoplot(term_aln, xlab="Reference genome position",
			geom = "rect", aes(fill = strand, colour = strand), show.legend = FALSE) +
 			facet_grid(strand ~ ., labeller = labeller(strand = strand.labs)) +
  			theme_calc() + 
			geom_vline(xintercept=plus_term_circ3, linetype="dashed", colour="springgreen3", linewidth=1) + 
  			geom_vline(xintercept=minus_term_circ3, linetype="dashed", colour="violet", linewidth=1)
		png("ggterm.png", width = 6, height = 8, units = "in", res = 300)
		print(ggterm)
		dev.off()
	}

	ggtau <- ggplot(tau_stats) +
		theme_calc() + 
		geom_point(data=subset(tau_stats, strand == "f"), aes(x=pos_adj, y=avg_tau, colour="plus")) +
		geom_point(data=subset(tau_stats, strand == "r"), aes(x=pos_adj, y=avg_tau, colour="minus")) +
		geom_label_repel(data=subset(tau_stats, strand == "f" & pos_adj == plus_term), 
                   aes(x=pos_adj, y=avg_tau, label=pos_adj), colour="#1B9E77",
                   show.legend = FALSE) + 
		geom_label_repel(data=subset(tau_stats, strand == "r" & pos_adj == minus_term), 
                   aes(x=pos_adj, y=avg_tau, label=pos_adj), colour="#7570B3",
                   show.legend = FALSE) +
		labs(x = "Reference genome position",
			y = "tau",
			colour = "strand") +
		scale_color_manual(values = colours) +
		scale_x_continuous(labels = comma) +
		scale_y_continuous(limits=c(0,1))

	ggdepth <- ggplot() +
 		geom_vline(xintercept=plus_term, linetype="dashed", colour="springgreen3", linewidth=1) +
  		geom_vline(xintercept=minus_term, linetype="dashed", colour="violet", linewidth=1) +
  		geom_line(data=sum_cov, aes(colour="both", x=pos, y=rollmean(covs, window, na.pad = TRUE, align = "right"))) +
  		geom_line(data=subset(tau, strand == "f"), aes(colour="plus", x=pos, y=rollmean(cov, window, na.pad = TRUE, align = "right"))) +
  		geom_line(data=subset(tau, strand == "r"), aes(colour="minus", x=pos, y=rollmean(cov, window, na.pad = TRUE, align = "right"))) +
  		theme_calc() +
  		labs(x = "Reference genome position",
       		y = "Read depth",
      		colour = "Strand") +
  		scale_color_manual(values = colours) +
  		scale_x_continuous(labels = comma) +
  		guides(colour = guide_legend(override.aes = list(linewidth = 3)))

	report <- read_docx() %>%
		body_add_par(value = paste("NanoTerm Report: ", name), style = "heading 1") %>%
		body_add_par("", style = "Normal") %>%
		body_add_par("Run details", style = "heading 2") %>%
		body_add_par(value = paste("Generated on: ", date, sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Reference genome: ", "$refName", sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Reference genome length: ", len, sep = ""), style = "Normal") %>%
		body_add_par("Alignment details", style = "heading 2") %>%
		body_add_par(value = paste("Sequencing platform: ", "$params.seqplat", sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Number of sequence reads: ", totalReads, sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Number of reads mapped: ", mappedReads, " (", percentMapped, "%)", sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Number of reads not mapped: ", unmappedReads, " (", percentUnmapped, "%)", sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Average read length: ", aveReadLen, sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Maximum read length: ", maxReadLen, sep = ""), style = "Normal") %>%
		body_add_par("Phage prediction", style = "heading 2") %>%
		body_add_par(value = paste("The predicted termini are ", location, " within the reference genome.", sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Plus terminus: ", plus_term, sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Minus terminus: ", minus_term, sep = ""), style = "Normal") %>%
		body_add_par(value = paste("The distance between predicted termini is ", term_dist, " nucleotides.", sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Class: ", class, sep = ""), style = "Normal") %>%
		body_add_par(value = paste("Subclass: ", subclass, sep = ""), style = "Normal") %>%
		body_add_par("", style = "Normal") %>%
		body_add_table(table, style = "table_template", first_column = TRUE) %>%
		body_add_par("", style = "Normal") %>%
		body_add_break(pos = "after") %>%
		body_add_par("Figures", style = "heading 2") %>%
		body_add_par("", style = "Normal") %>%
		body_add_gg(value = ggtau, style = "centered", height = 3.25) %>%
		body_add_par(value = "Figure 1. The tau value calculated for each genome position.") %>%
		body_add_par("", style = "Normal") %>%
		body_add_par("", style = "Normal") %>%
		body_add_gg(value = ggdepth, style = "centered", height = 3.25) %>%
		body_add_par(value = "Figure 2. The total read depth of the sequencing run, graphed as a rolling average with a window size equal to 1% of the reference genome length.") %>%
 		body_add_par("", style = "Normal")
	if (terminalReads == TRUE & location == "internal"){
		report <- body_add_img(report, src = "ggterm.png", style = "centered", width = 6, height = 8) %>%
		body_add_par(value = "Figure 3. Reads that cover part or all of the region between the predicted termini.")
	}
	if (terminalReads == TRUE & location == "terminal"){
		report <- body_add_img(report, src = "ggterm.png", style = "centered", width = 6, height = 8) %>%
		body_add_par(value = "Figure 3. Reads that cover part or all of the region between the predicted termini.  Due to the terminal nature of the predicted termini, the read are mapped against the third circular permutation of the reference genome.")
	}
	print(report, target = "./report.docx")
	"""
}

// This process converts to the .docx report to a PDF
process doc2pdf {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true
	
	input:
		path report
		
	output:
		path 'report.pdf'
	
	script:
	"""
	pandoc -V geometry:margin=1in -f docx -t latex -o report.pdf $report 
	"""	
}

// This process extracts the DTR or COS sequences
// It also rearranges the genome based on the predicted termini
process fastaOut {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true

	input:
		path classification
		path refseq
		path circular_permutation1
		path circular_permutation2
		path circular_permutation3
		path circular_permutation4
		path circular_permutation5

	output:
		path 'rearranged_genome.fasta'
		path 'DTR_sequence.fasta', optional: true

	script:
		"""
		#!/usr/bin/env python3

		import pandas
		df = pandas.read_csv('classification.csv')
		plus_term = df.iat[0,5]
		minus_term = df.iat[0,7]
		type = df.iat[0,11]
		subclass = df.iat[0,10]
		location = df.iat[0,9]
		plus_term_circ3 = df.iat[0,13]
		minus_term_circ3 = df.iat[0,14]

		if location == "internal" and type == "DTR":
			f = open('refseq.fasta', 'r')
			next(f)
			for line in f:
				refseq = str(line.rstrip())
			f.close()

			newseq = refseq[plus_term:] + refseq[:minus_term]

			sourceFile = open('rearranged_genome.fasta', 'w')
			print('>rearranged_reference_genome', file = sourceFile)
			print(newseq, file = sourceFile)
			sourceFile.close()
		
			DTRseq = refseq[plus_term:minus_term]

			sourceFile = open('DTR_sequence.fasta', 'w')
			print('>DTR_sequence', file = sourceFile)
			print(DTRseq, file = sourceFile)
			sourceFile.close()

		if location == "terminal" and type == "DTR":
			f = open('circular_permutation3.fasta', 'r')
			next(f)
			for line in f:
				refseq = str(line.rstrip())
			f.close()

			newseq = refseq[plus_term_circ3:] + refseq[:minus_term_circ3]

			sourceFile = open('rearranged_genome.fasta', 'w')
			print('>rearranged_reference_genome', file = sourceFile)
			print(newseq, file = sourceFile)
			sourceFile.close()
						
			DTRseq = refseq[plus_term_circ3:minus_term_circ3]

			sourceFile = open('DTR_sequence.fasta', 'w')
			print('>DTR_sequence', file = sourceFile)
			print(DTRseq, file = sourceFile)
			sourceFile.close()

		if location == "internal" and type == "COS":
			f = open('refseq.fasta', 'r')
			next(f)
			for line in f:
				refseq = str(line.rstrip())
			f.close()

			newseq = refseq[plus_term:] + refseq[:minus_term]

			sourceFile = open('rearranged_genome.fasta', 'w')
			print('>rearranged_reference_genome', file = sourceFile)
			print(newseq, file = sourceFile)
			sourceFile.close()
		
		if location == "terminal" and type == "COS":
			f = open('circular_permutation3.fasta', 'r')
			next(f)
			for line in f:
				refseq = str(line.rstrip())
			f.close()

			newseq = refseq[plus_term_circ3:] + refseq[:minus_term_circ3]

			sourceFile = open('rearranged_genome.fasta', 'w')
			print('>rearranged_reference_genome', file = sourceFile)
			print(newseq, file = sourceFile)
			sourceFile.close()
		"""
}

workflow {
	allseq_ch = catFastq(seq_ch)
	raw_ch = rawSeq(ref_ch)
	len_ch = refLen(raw_ch)
	refseq_ch = refSeq(ref_ch)
	circ_ch = permute(refseq_ch)
	aln_ch = mapping(plat_ch, refseq_ch, allseq_ch, circ_ch)
	alnstats_ch = alignStats(aln_ch, ref_ch)
	sep_ch = strandSep(aln_ch)
	bed_ch = bed(sep_ch)
	cov_ch = cov(sep_ch)
	tau_ch = tau(len_ch, bed_ch, cov_ch)
	class_ch = classify(len_ch, tau_ch, name_ch, alnstats_ch, seq_ch, ref_ch, plat_ch)
	termreads_ch = terminalReads(aln_ch, class_ch)
	report_ch = report(ref_ch, name_ch, plat_ch, alnstats_ch, class_ch, tau_ch, termreads_ch)
	doc2pdf(report_ch)
	fastaOut(refseq_ch, class_ch, circ_ch)
}