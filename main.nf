#!/usr/bin/env nextflow

ref_ch = Channel.of(params.ref)
seq_ch = Channel.of(params.seq)

println "\nReference: $params.ref\nSequence reads: $params.seq\n"

process rawSeq {
	input:
		path ref
	
	output:
		path 'fasta.txt'
		
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
		path 'extref.fasta'	
	
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

process extRawSeq {
	input:
		path extref
	
	output:
		path 'extFasta.txt'
		
	script:
	"""
	grep -v ">" $extref | tr -d "\\n" > extFasta.txt
	"""
}

process mapping {
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

process strandSep {
	input:
		path aln
		
	output:
		path 'aln_f_sorted.bam'
		path 'aln_r_sorted.bam'

	script:
		"""
		samtools view -F 16 -b $aln > aln_f.bam
		samtools sort aln_f.bam > aln_f_sorted.bam
		
		samtools view -f 16 -b $aln > aln_r.bam
		samtools sort aln_r.bam > aln_r_sorted.bam
		"""
}

process bed {
	input:
		path aln_f_sorted
		path aln_r_sorted
		
	output:
		path 'aln_f.bed'
		path 'aln_r.bed'
		
	script:
	"""
	bamToBed -i $aln_f_sorted > aln_f.bed
	bamToBed -i $aln_r_sorted > aln_r.bed
	"""
}

process cov {	
	input:
		path aln_f_sorted
		path aln_r_sorted
		
	output:
		path 'aln_f_cov.txt'
		path 'aln_r_cov.txt'
	
	script:
	"""
	samtools mpileup -d 0 -o aln_f_cov.txt $aln_f_sorted
	samtools mpileup -d 0 -o aln_r_cov.txt $aln_r_sorted
	"""
}

process tau {
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

process extMapping {
	input:
		path seq
		path extref

	output:
		path 'ext_aln.sam'
		
	script:
		"""
		minimap2 -ax map-ont $extref $params.seq > ext_aln.sam
		"""
}

process extStrandSep {
	input:
		path ext_aln
		
	output:
		path 'ext_aln_f_sorted.bam'
		path 'ext_aln_r_sorted.bam'

	script:
		"""
		samtools view -F 16 -b $ext_aln > ext_aln_f.bam
		samtools sort ext_aln_f.bam > ext_aln_f_sorted.bam
		
		samtools view -f 16 -b $ext_aln > ext_aln_r.bam
		samtools sort ext_aln_r.bam > ext_aln_r_sorted.bam
		"""
}

process extBed {
	input:
		path ext_aln_f_sorted
		path ext_aln_r_sorted
		
	output:
		path 'ext_aln_f.bed'
		path 'ext_aln_r.bed'
		
	script:
	"""
	bamToBed -i $ext_aln_f_sorted > ext_aln_f.bed
	bamToBed -i $ext_aln_r_sorted > ext_aln_r.bed
	"""
}

process extCov {	
	input:
		path ext_aln_f_sorted
		path ext_aln_r_sorted
		
	output:
		path 'ext_aln_f_cov.txt'
		path 'ext_aln_r_cov.txt'
	
	script:
	"""
	samtools mpileup -d 0 -o ext_aln_f_cov.txt $ext_aln_f_sorted
	samtools mpileup -d 0 -o ext_aln_r_cov.txt $ext_aln_r_sorted
	"""
}

process extTau {
	input:
		path ext_aln_f
		path ext_aln_r
		path ext_aln_f_cov
		path ext_aln_r_cov
		path extFasta
		
	output:
		path 'ext_tau.csv'
		
	script:
	"""
	#!/usr/bin/env Rscript

	f_bed <- read.table("$ext_aln_f", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	f_bed <- setNames(f_bed, c("chr","start","end"))
	f_bed['start'] <- f_bed['start'] + 1

	r_bed <- read.table("$ext_aln_r", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric","NULL","NULL","NULL"))
	r_bed <- setNames(r_bed, c("chr","start","end"))
	r_bed['start'] <- r_bed['start'] + 1

	f_cov <- read.table("$ext_aln_f_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","NULL","numeric","NULL","NULL"))
	f_cov <- setNames(f_cov, c("chr","pos","cov"))
	f_cov['strand'] <- "f"

	r_cov <- read.table("$ext_aln_r_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","NULL","numeric","NULL","NULL"))
	r_cov <- setNames(r_cov, c("chr","pos","cov"))
	r_cov['strand'] <- "r"
	
	fa <- scan(file="$extFasta", what="string")
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
	
	ext_tau <- rbind(f_tau, r_tau)

	ext_tau['tau'] <- ext_tau['SPC'] / ext_tau['cov']
	
	write.csv(ext_tau, "ext_tau.csv")
	"""	
}

workflow {
	len_ch = rawSeq(ref_ch)
	ext_ch = extendRef(len_ch, ref_ch)
	aln_ch = mapping(seq_ch, ref_ch)
	sep_ch = strandSep(aln_ch)
	bed_ch = bed(sep_ch)
	cov_ch = cov(sep_ch)
	tau(bed_ch, cov_ch, len_ch)
	extlen_ch = extRawSeq(ext_ch)
	extaln_ch = extMapping(seq_ch, ext_ch)
	extsep_ch = extStrandSep(extaln_ch)
	extbed_ch = extBed(extsep_ch)
	extcov_ch = extCov(extsep_ch)
	extTau(extbed_ch, extcov_ch, extlen_ch)
}