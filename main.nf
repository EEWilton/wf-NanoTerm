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
	grep -v ">" $reference | tr -d "\\n" > rawseq.txt
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
	grep ">" $reference > refseq.fasta
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
		path all
		path refseq
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
			minimap2 -ax map-ont $params.reference $all > aln.sam
			minimap2 -ax map-ont $circular_permutation1 $all  > aln_circ1.sam
			minimap2 -ax map-ont $circular_permutation2 $all  > aln_circ2.sam
			minimap2 -ax map-ont $circular_permutation3 $all  > aln_circ3.sam
			minimap2 -ax map-ont $circular_permutation4 $all > aln_circ4.sam
			minimap2 -ax map-ont $circular_permutation5 $all  > aln_circ5.sam
			"""
		
		else if( seqplat == 'illumina' )
			"""
			minimap2 -asr $params.reference $params.fastq > aln.sam
			minimap2 -asr $circular_permutation1 $all > aln_circ1.sam
			minimap2 -asr $circular_permutation2 $all  > aln_circ2.sam
			minimap2 -asr $circular_permutation3 $all > aln_circ3.sam
			minimap2 -asr $circular_permutation4 $all  > aln_circ4.sam
			minimap2 -asr $circular_permutation5 $pall  > aln_circ5.sam
			"""

		else 
			 error "Invalid sequencing platform: $params.seqplat"
}

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
		
	script:
	"""
	#!/usr/bin/env Rscript

	len <- $genLen
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
	write.csv(tau_circ5, "tau_circ5.csv")
	
	"""	
}

process stats {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true

	input:
		path tau
		path tau_circ1
		path tau_circ2
		path tau_circ3
		path tau_circ4
		path tau_circ5
		
	output:
		path 'stats.csv'
	
	script:
	"""
	#! /usr/bin/env python3

	import heapq
	import itertools
	import pandas as pd
	import numpy as np
	from sklearn.tree import DecisionTreeRegressor
	from scipy import stats
	from statsmodels.sandbox.stats.multicomp import multipletests

	def removePeaks(arr,n):
		arr=np.array(arr)
		peak_pos=arr.argsort()[-n:][::-1]
		arr2=np.delete(arr,peak_pos)
		return arr2

	def gamma(X):
		X = np.array(X, dtype=np.int64)
		v = removePeaks(X, 3)

		dist_max = float(max(v))
		if dist_max == 0:
			return np.array([1.00] * len(X))

		actual = np.bincount(v)
		fit_alpha, fit_loc, fit_beta = stats.gamma.fit(v)
		expected = stats.gamma.pdf(np.arange(0, dist_max + 1, 1), fit_alpha, loc=fit_loc, scale=fit_beta) * sum(actual)

		return stats.gamma.pdf(X, fit_alpha, loc=fit_loc, scale=fit_beta)

	def peaksDecisionTree(read_depth, starting_pos_depth, tau):
		L = len(read_depth[0])
		res = pd.DataFrame({"Position": np.array(range(L)) + 1, "SPC_plus": starting_pos_depth[0],
                        "tau_plus": tau[0], "tau_minus": tau[1], "SPC_minus": starting_pos_depth[1],
                        "cov_plus": read_depth[0], "cov_minus": read_depth[1]})

		res["cov"] = res["cov_plus"].values + res["cov_minus"].values

		res["R_plus"] = list(map(float, starting_pos_depth[0])) // np.mean(starting_pos_depth[0])
		res["R_minus"] = list(map(float, starting_pos_depth[1])) // np.mean(starting_pos_depth[1])

		regr = DecisionTreeRegressor(max_depth=3, min_samples_leaf=100)
		X = np.arange(L)
		X = X[:, np.newaxis]
		y = res["cov"].values
		regr.fit(X, y)

		y_1 = regr.predict(X)
		res["covnode"] = y_1
		covnodes = np.unique(y_1)
		thres = np.mean(read_depth[0]) / 2
		covnodes = [n for n in covnodes if n > thres]

		for node in covnodes:
			X = res[res["covnode"] == node]["SPC_plus"].values
			res.loc[res["covnode"] == node, "pval_plus"] = gamma(X)
			X = res[res["covnode"] == node]["SPC_minus"].values
			res.loc[res["covnode"] == node, "pval_minus"] = gamma(X)

		res.loc[res.pval_plus > 1, 'pval_plus'] = 1.00
		res.loc[res.pval_minus > 1, 'pval_minus'] = 1.00
		res = res.fillna(1.00)

		res['pval_plus_adj'] = multipletests(res["pval_plus"].values, alpha=0.01, method="bonferroni")[1]
		res['pval_minus_adj'] = multipletests(res["pval_minus"].values, alpha=0.01, method="bonferroni")[1]

		res = res.fillna(1.00)

		res_plus = pd.DataFrame(
			{"Position": res['Position'], "tau": res['tau_plus'], "pval_gamma": res['pval_plus'], "pval_gamma_adj": res['pval_plus_adj']})
		res_minus = pd.DataFrame(
			{"Position": res['Position'], "tau": res['tau_minus'], "pval_gamma": res['pval_minus'], "pval_gamma_adj": res['pval_minus_adj']})

		res_plus.sort_values("tau", ascending=False, inplace=True)
		res_minus.sort_values("tau", ascending=False, inplace=True)
	
		res_plus.reset_index(drop=True, inplace=True)
		res_minus.reset_index(drop=True, inplace=True)

		return res, res_plus, res_minus

	def selectSignificant(table, pvalue, limit):
		table_pvalue = table.loc[lambda df: df.pval_gamma_adj < pvalue, :]
		table_pvalue_limit = table_pvalue.loc[lambda df: df.tau < limit, :]
		table_pvalue_limit.reset_index(drop=True, inplace=True)
		return table_pvalue_limit
	
	data = pd.read_csv("$tau")

	f = data[data['strand'] == "f"]
	r = data[data['strand'] == "r"]

	f_cov = f['cov'].to_numpy()
	r_cov = r['cov'].to_numpy()

	read_depth = np.array([f_cov, r_cov])

	f_spc = f['SPC'].to_numpy()
	r_spc = r['SPC'].to_numpy()

	starting_pos_depth = np.array([f_spc, r_spc])

	f_tau = f['tau'].to_numpy()
	r_tau = r['tau'].to_numpy()

	tau = np.array([f_tau, r_tau])
	
	surrounding = 20
	gen_len = len(f_cov)

	phage_norm, phage_plus_norm, phage_minus_norm = peaksDecisionTree(read_depth, starting_pos_depth, tau)
	
	plus_significant = selectSignificant(phage_plus_norm, 1.0 / gen_len, 1.0)
	minus_significant = selectSignificant(phage_minus_norm, 1.0 / gen_len, 1.0)

	phage_norm.to_csv("stats.csv")
	"""
}

process report {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true
	
	input:
		path fastq
		path reference
		path stats
		val name
		val seqplat
		val totalReads
		val mappedReads
		val unmappedReads
		val aveReadLen
		val maxReadLen
		
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

	data <- read.csv("stats.csv", quote="\\"", comment.char="")

	meanDepth <- summarise(data, meanDepth = mean(cov))
	
	print(meanDepth)
	
	date <- format(Sys.Date())
	name <- as.character("$name")
	totalReads <- as.integer("$totalReads")
	mappedReads <- as.integer("$mappedReads")
	unmappedReads <- as.integer("$unmappedReads")
	aveReadLen <- as.integer("$aveReadLen")
	maxReadLen <- as.integer("$maxReadLen")
	
	len <- max(data['Position'])
	window <- 0.01 * len
   
	max_plus <- max(data['cov_plus'])
	max_minus <- max(data['cov_minus'])

	top_tau_plus <- data %>%
		filter(tau_plus < 1) %>%
		filter(tau_plus > 0.1) %>%
		arrange(desc(tau_plus)) %>%
		filter(pval_plus_adj <= 0.5) %>%
		filter(cov_plus > (0.05 * max_plus)) 
	
	top_tau_minus <- data %>%
		filter(tau_minus < 1) %>%
		filter(tau_minus > 0.1) %>%
		arrange(desc(tau_minus)) %>%
		filter(pval_minus_adj <= 0.5) %>%
		filter(cov_minus > (0.05 * max_minus)) 

	f_term <- top_tau_plus[1,2]
	r_term <- top_tau_minus[1,2]
	
	if (is.na(f_term)){
		plus_term <- "NA"
	} else {
		plus_term <- f_term
	}

	if (is.na(r_term)){
		minus_term <- "NA"
	} else {
		minus_term <- r_term
	}

	within_DTR <- data %>%
		filter(Position > plus_term) %>%
		filter(Position < minus_term)

	outside_DTR <- data %>%
		filter(Position < plus_term | Position > minus_term)

	mean_DTR_depth <-  summarise(within_DTR, mean_DTR_depth = mean(cov))
	mean_notDTR_depth <- summarise(outside_DTR, mean_notDTR_depth = mean(cov))
	
	DTR_depth_ratio <- mean_DTR_depth / mean_notDTR_depth

	if (is.integer(plus_term)){
		f_term_tau <- top_tau_plus[1,4]
		f_term_p <- top_tau_plus[1,15]
	} else {
		f_term_tau <- "NA"
		f_term_p <- "NA"
	}

	if (is.integer(minus_term)){
		r_term_tau <- top_tau_minus[1,5]
		r_term_p <- top_tau_minus[1,16]
	} else {
		r_term_tau <- "NA"
		r_term_p <- "NA"
	}

	Strand <- c("+","-")
	Terminus <- c(plus_term, minus_term)
	tau <- c(f_term_tau, r_term_tau)
	p_value <- c(f_term_p, r_term_p)
	table <- data.frame(Strand, Terminus, tau, p_value)

	
	if (is.integer(plus_term) && is.integer(minus_term)) {
		term_dist <- minus_term - plus_term
	} else {
		term_dist = "NA"
	}
	
	colours <- c("plus" = "springgreen4", "minus" = "purple")

	if (is.integer(term_dist) == TRUE){
		if (term_dist > 20){
			class <- "DTR"
		} else if (term_dist < 0){
			class <- "COS 3′"
		} else {
			class <- "COS 5′"
	}} else {
		if (is.integer(plus_term) | is.integer(minus_term)) {
			class <- "Headful"
		} else {
			class <- "Headful with no pac site, mu-like, or other"
	}}

	if (class == "Headful"){
		if (is.integer(plus_term)){
			pack_dir <- "Forward"
		} else if (is.integer(minus_term)){
			pack_dir <- "Reverse"
		}
	}
	
	if (class == "Headful"){
		if (is.integer(plus_term)){
			concat <- (1 - f_term_tau)/f_term_tau
		} else if (is.integer(minus_term)){
			concat <- (1 - r_term_tau)/r_term_tau
		}
	}
	
	depth <- ggplot(data = data, aes(x=Position, y=rollmean(cov, window, na.pad = TRUE, align = "right")))
		if (is.integer(plus_term)) depth <- depth + geom_vline(xintercept=plus_term, linetype="dashed", colour="springgreen3", linewidth=1.1) 
		if (is.integer(minus_term)) depth <- depth + geom_vline(xintercept=minus_term, linetype="dashed", colour="violet", linewidth=1.1) 
	depth <- depth + geom_line(data=data, aes(x=Position, y=rollmean(cov_plus, window, na.pad = TRUE, align = "right"), colour="plus")) +
		geom_line(data=data, aes(x=Position, y=rollmean(cov_minus, window, na.pad = TRUE, align = "right"), colour="minus")) +
		geom_line() +
		labs(x = "Reference genome position",
			y = "Read depth",
			colour = "Strand") +
		scale_color_manual(values = colours) +
		scale_x_continuous(labels = comma) +
		scale_y_continuous(limits = c(0, NA)) +
		guides(colour = guide_legend(override.aes = list(linewidth = 3))) +
		theme_calc()

	tau <- ggplot(data=data) +
		theme_calc() + 
		geom_point(data=subset(data, tau_plus < 1), aes(x=Position, y=tau_plus, colour="plus")) +
		geom_point(data=subset(data, tau_minus < 1), aes(x=Position, y=tau_minus, colour="minus")) +
		geom_label_repel(data=subset(data, Position == plus_term), 
                   aes(x=Position, y=tau_plus,label=Position), colour="springgreen4",
                   show.legend = FALSE) + 
		geom_label_repel(data=subset(data,  Position == minus_term), 
                   aes(x=Position, y=tau_minus,label=Position), colour="purple",
                   show.legend = FALSE) +
		labs(x = "Reference genome position",
			y = "tau",
			colour = "Strand") +
		scale_color_manual(values = colours) +
		scale_x_continuous(labels = comma) +
		scale_y_continuous(limits=c(0,1))
	
	report <- read_docx() %>%
		body_add_par(value = paste("NanoTerm Report: ", name), style = "heading 1") %>%
		body_add_par("", style = "Normal") %>%
		body_add_par("Run details", style = "heading 2") %>%
		body_add_par(value = paste("Generated on: ", date), style = "Normal") %>%
		body_add_par(value = paste("Input sequences: ", "$params.fastq"), style = "Normal") %>%
		body_add_par(value = paste("Reference genome: ", "$params.reference"), style = "Normal") %>%
		body_add_par(value = paste("Reference genome length: ", len), style = "Normal") %>%
		body_add_par("Alignment details", style = "heading 2") %>%
		body_add_par(value = paste("Sequencing platform:", "$params.seqplat"), style = "Normal") %>%
		body_add_par(value = paste("Number of sequence reads: ", totalReads), style = "Normal") %>%
		body_add_par(value = paste("Number of reads aligned: ", mappedReads), style = "Normal") %>%
		body_add_par(value = paste("Number of reads not aligned: ", unmappedReads), style = "Normal") %>%
		body_add_par(value = paste("Average read length: ", aveReadLen), style = "Normal") %>%
		body_add_par(value = paste("Maximum read length: ", maxReadLen), style = "Normal") %>%
		body_add_par(value = paste("Average read depth: ", as.integer(meanDepth)), style = "Normal") %>%
		body_add_par("", style = "Normal") %>%
		body_add_par("Phage prediction", style = "heading 2") %>%
		body_add_table(table, style = "table_template", first_column = TRUE) %>%
		body_add_par("", style = "Normal") %>%
		body_add_par(value = paste("Phage class: ", class), style = "Normal")

	if (class == "COS 5′"){
		report <- body_add_par(report, value = paste("Cohesive sequence: ", class), style = "Normal")
	} else if (class == "COS 3′"){
		report <- body_add_par(report, value = paste("Cohesive sequence: ", class), style = "Normal")
	} else if (class == "DTR"){
		report <- body_add_par(report, value = paste("DTR length: ", term_dist), style = "Normal")
		report <- body_add_par(report, value = paste("Average read depth within DTR: ", mean_DTR_depth), style = "Normal")
		report <- body_add_par(report, value = paste("DTR read depth / non-DTR read depth: ", DTR_depth_ratio), style = "Normal")
		report <- body_add_par(report, value = "If there is more the 2x the read depth within the predicted DTR, this is consistent with DTR packaging", style = "Normal")
	} else if (class == "headful"){
		report <- body_add_par(report, value = paste("Packaging direction: ", pack_dir), style = "Normal")
		report <- body_add_par(report, value = paste("Number of genome copies per concatamer: ", concat), style = "Normal")
	} else {
	}

	report <- body_add_par(report, "Figures", style = "heading 2") %>%
		body_add_par("", style = "Normal") %>%
		body_add_gg(value = tau, style = "centered", height = 3.25) %>%
		body_add_par(value = "Figure 1. The tau value calculated for each genome position.") %>%
		body_add_par("", style = "Normal") %>%
		body_add_gg(value = depth, style = "centered", height = 3.25) %>%
		body_add_par(value = "Figure 2. The total read depth of the sequencing run, graphed as a rolling average with a window size equal to 1% of the reference genome length.  Black is the sum of forward and reverse read depth.")
 
	print(report, target = "./report.docx")
	"""
}

process doc2pdf {
	publishDir "${params.out_dir}", mode: 'copy', overwrite: true
	
	input:
		path report
		path out_dir
		
	output:
		path 'report.pdf', optional: true
	
	script:
	"""
	libreoffice --headless --convert-to pdf  $report --outdir $params.out_dir
	"""	
}

workflow {
	allseq_ch = catFastq(seq_ch)
	raw_ch = rawSeq(ref_ch)
	len_ch = refLen(raw_ch)
	refseq_ch = refSeq(ref_ch)
	circ_ch = permute(refseq_ch)
	aln_ch = mapping(plat_ch, allseq_ch, refseq_ch, circ_ch)
	alnstats_ch = alignStats(aln_ch)
	sep_ch = strandSep(aln_ch)
	bed_ch = bed(sep_ch)
	cov_ch = cov(sep_ch)
	tau_ch = tau(len_ch, bed_ch, cov_ch)
	stats_ch = stats(tau_ch)
	report_ch = report(stats_ch, ref_ch, refseq_ch, name_ch, plat_ch, alnstats_ch)
	doc2pdf(report_ch, outdir_ch)
}