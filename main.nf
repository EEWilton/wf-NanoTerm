#!/usr/bin/env nextflow

ref_ch = Channel.of(params.ref)
seq_ch = Channel.of(params.seq)

println "\nReference: $params.ref\nSequence reads: $params.seq\n"

process raw_seq {
	input:
		path ref
	
	output:
		path "fasta.txt"
		
	script:
	"""
	grep -v ">" $ref | tr -d "\\n" > fasta.txt
	"""
}

process extendRef {
	input:
		path fasta
		path ref
			
	output:
		path "extref.fasta"
	
	script:
	"""
	#!/usr/bin/env python3
	
	f = open('fasta.txt', 'r')
	refseq = f.read()
	f.close()
	
	edge = 500
	
	extref = refseq[-edge:] + refseq + refseq[:edge]

	sourceFile = open('extref.fasta', 'w')
	print('>extended_reference', file = sourceFile)
	print(extref, file = sourceFile)
	sourceFile.close()
	"""
}

process runMinimap2 {
	conda 'my-env.yaml'

	input:
		path seq
		path ref

	output:
		path 'aln.sam'
		
	script:
		"""
		minimap2 -ax map-ont $params.ref $params.seq > aln.sam
		"""
}

process Fstrand {
	conda './my-env.yaml'
	
	input:
		path aln
		
	output:
		path 'aln_f_sorted.bam'

	script:
		"""
		samtools view -F 16 -b $aln > aln_f.bam
		samtools sort aln_f.bam > aln_f_sorted.bam
		"""
}

process Rstrand {
	conda './my-env.yaml'
	
	input:
		path aln
		
	output:
		path 'aln_r_sorted.bam'

	script:
		"""
		samtools view -f 16 -b $aln > aln_r.bam
		samtools sort aln_r.bam > aln_r_sorted.bam
		"""
}

process Fbed {
	conda './my-env.yaml'
	
	input:
		path aln_f_sorted
		
	output:
		path 'aln_f.bed'
		
	script:
	"""
	bamToBed -i $aln_f_sorted > aln_f.bed
	"""
}

process Rbed {
	conda './my-env.yaml'
	
	input:
		path aln_r_sorted
		
	output:
		path 'aln_r.bed'
		
	script:
	"""
	bamToBed -i $aln_r_sorted > aln_r.bed
	"""
}

process Fcov {
	conda './my-env.yaml'
	
	input:
		path aln_f_sorted
		
	output:
		path 'aln_f_cov.txt'
	
	script:
	"""
	samtools mpileup -d 0 -o aln_f_cov.txt $aln_f_sorted
	"""
}

process Rcov {
	conda './my-env.yaml'
	
	input:
		path aln_r_sorted
		
	output:
		path 'aln_r_cov.txt'
	
	script:
	"""
	samtools mpileup -d 0 -o aln_r_cov.txt $aln_r_sorted
	"""
}

process tau_R {
	conda './my-env.yaml'
	input:
		path aln_f
		path aln_r
		path aln_f_cov
		path aln_r_cov
		path fasta
		
	output:
		path "tau.csv"
		
	script:
	"""
	#!/usr/bin/env Rscript

	f_bed <- read.table("$aln_f", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	f_bed <- setNames(f_bed, c("chr","start","end"))
	f_bed['start'] <- f_bed['start'] + 1

	r_bed <- read.table("$aln_r", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	r_bed <- setNames(r_bed, c("chr","start","end"))
	r_bed['start'] <- r_bed['start'] + 1

	f_cov <- read.table("$aln_f_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","NULL","numeric","NULL","NULL"))
	f_cov <- setNames(f_cov, c("chr","pos","cov"))
	f_cov['strand'] <- "f"

	r_cov <- read.table("$aln_r_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","NULL","numeric","NULL","NULL"))
	r_cov <- setNames(r_cov, c("chr","pos","cov"))
	r_cov['strand'] <- "r"
	
	fa <- scan(file="$fasta", what="string")
	len <- nchar(fa)
	pos <- seq(from = 1,to = len, by = 1)
	
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
	"""	
}


workflow {
	len_ch = raw_seq(ref_ch)
	ext_ch = extendRef(len_ch, ref_ch)
	aln_ch = runMinimap2(seq_ch, ref_ch)
	F_ch = Fstrand(aln_ch)
	R_ch = Rstrand(aln_ch)
	Fbed_ch = Fbed(F_ch)
	Rbed_ch = Rbed(R_ch)
	Fcov_ch = Fcov(F_ch)
	Rcov_ch = Rcov(R_ch)
	tau_R(Fbed_ch, Rbed_ch, Fcov_ch, Rcov_ch, len_ch)
}