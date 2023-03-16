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
		path extref

	output:
		path 'aln.sam'
		
	script:
		"""
		minimap2 -ax map-ont $extref $params.seq > aln.sam
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

process bedDF {
	conda './my-env.yaml'

	input:
		path aln_f
		path aln_r
	
	output:
		path "f_bed.txt"
		path "r_bed.txt"
		
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

	sink("f_bed.txt")
	print(f_bed)
	sink()
	
	sink("r_bed.txt")
	print(r_bed)
	sink()
	"""
}


process covDF {
	conda './my-env.yaml'

	input:
		path aln_f_cov
		path aln_r_cov
	
	output:
		path "f_cov.txt"
		path "r_cov.txt"
		
	script:
	"""
	#!/usr/bin/env Rscript
	
	f_cov <- read.table("$aln_f_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","NULL","numeric","NULL","NULL"))
	f_cov <- setNames(f_cov, c("chr","pos","cov"))
	f_cov['strand'] <- "f"

	r_cov <- read.table("$aln_r_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","NULL","numeric","NULL","NULL"))
	r_cov <- setNames(r_cov, c("chr","pos","cov"))
	r_cov['strand'] <- "r"
	
	sink("f_cov.txt")
	print(f_cov)
	sink()
	
	sink("r_cov.txt")
	print(r_cov)
	sink()
	"""
}

process spcDF {
	conda './my-env.yaml'

	input:
		path f_bed
		path r_bed
		path extref
	
	output:
		path "f_spc.txt"
		path "r_spc.txt"
		
	script:
	"""
	#!/usr/bin/env Rscript
	
	fa <- scan(file="$extref", what="string")
	len <- nchar(fa)
	pos <- seq(1,len,1)
	
	f_spc <- data.frame(pos = pos, SPC = 'NA')
	f_bed <- read.csv("$f_bed", sep="")

	for (nt in pos){
		f_spc[nt,'SPC'] <- length(which(f_bed['start']==nt))
		}

	sink("f_spc.txt")
	print(f_spc)
	sink()
	
	r_spc <- data.frame(pos = pos, SPC = 'NA')
	r_bed <- read.csv("$r_bed", sep="")

	for (nt in pos){
		r_spc[nt,'SPC'] <- length(which(r_bed['end']==nt))
		}

	sink("r_spc.txt")
	print(r_spc)
	sink()
	"""
}
process tauDF {
	conda './my-env.yaml'

	input:
	path f_cov
	path r_cov
	path f_spc
	path r_spc
	
	output:
	path "tau.txt"
	
	script:
	"""
	#!/usr/bin/env Rscript
	
	f_cov <- read.csv("$f_cov", sep="")
	f_spc <- read.csv("$f_spc", sep="")
	
	r_cov <- read.csv("$r_cov", sep="")
	r_spc <- read.csv("$r_spc", sep="")
	
	f_tau <- merge(f_cov, f_spc, by="pos")
	r_tau <- merge(r_cov, r_spc, by="pos")
	
	tau <- rbind(f_tau, r_tau)
	
	tau['tau'] <- tau['SPC']/tau['cov']
	
	sink("tau.txt")
	print(tau)
	sink()
	"""
}

workflow {
	len_ch = raw_seq(ref_ch)
	ext_ch = extendRef(len_ch, ref_ch)
	aln_ch = runMinimap2(seq_ch, ext_ch)
	F_ch = Fstrand(aln_ch)
	R_ch = Rstrand(aln_ch)
	Fbed_ch = Fbed(F_ch)
	Rbed_ch = Rbed(R_ch)
	Fcov_ch = Fcov(F_ch)
	Rcov_ch = Rcov(R_ch)
	bed_ch = bedDF(Fbed_ch, Rbed_ch)
	cov_ch = covDF(Fcov_ch, Rcov_ch)
	spc_ch = spcDF(bed_ch, len_ch)
	tau_ch = tauDF(cov_ch, spc_ch)
}