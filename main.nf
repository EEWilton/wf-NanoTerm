#!/usr/bin/env nextflow

ref_ch = Channel.of(params.fasta)
seq_ch = Channel.of(params.input)
plat_ch = Channel.of(params.seqplat)

println "\nReference: $params.fasta\nSequence reads: $params.input\n"

process rawSeq {
	input:
		path fasta
	
	output:
		path 'rawseq.txt'
		
	script:
	"""
	grep -v ">" $fasta | tr -d "\\n" > rawseq.txt
	"""
}

process extendRef {
	input:
		path rawseq
			
	output:
		path 'extref.fasta'	
	
	script:
	"""
	#!/usr/bin/env python3
	
	f = open('rawseq.txt', 'r')
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
		path input
		path fasta
		val seqplat

	output:
		path 'aln.sam'
		
	script:
		if( seqplat == 'nanopore' )
			"""
			minimap2 -ax map-ont $params.fasta $params.input > aln.sam
			"""
		
		else if( seqplat == 'illumina' )
			"""
			minimap2 -asr $params.fasta $params.input > aln.sam
			"""

		else 
			 error "Invalid sequencing platform: $params.seqplat"
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
	samtools depth -a -d 0 -o aln_f_cov.txt $aln_f_sorted
	samtools depth -a -d 0 -o aln_r_cov.txt $aln_r_sorted
	"""
}


process tau {
	publishDir "${params.outdir}", mode: 'copy', overwrite: true
	
	input:
		path aln_f
		path aln_r
		path aln_f_cov
		path aln_r_cov
		path fasta
		
	output:
		path 'tau.csv'
		
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
                    colClasses=c("character","numeric","numeric"))
	f_cov <- setNames(f_cov, c("chr","pos","cov"))
	f_cov['strand'] <- "f"

	r_cov <- read.table("$aln_r_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
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

process stats {
	publishDir "${params.outdir}", mode: 'copy', overwrite: true

	input:
		path tau
		
	output:
		path 'phage_norm.csv'
	
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

	def picMax(coverage, nbr_pic):
		picMaxPlus = heapq.nlargest(nbr_pic, zip(coverage[0], itertools.count()))
		picMaxMinus = heapq.nlargest(nbr_pic, zip(coverage[1], itertools.count()))
		TopFreqH = max(max(np.array(list(zip(*picMaxPlus))[0])), max(np.array(list(zip(*picMaxMinus))[0])))
		return picMaxPlus, picMaxMinus, TopFreqH

	def RemoveClosePicMax(picMax, gen_len, nbr_base):
		if nbr_base == 0:
			return picMax[1:], [picMax[0]]
		picMaxRC = picMax[:]
		posMax = picMaxRC[0][1]
		LimSup = posMax + nbr_base
		LimInf = posMax - nbr_base
		if LimSup < gen_len and LimInf >= 0:
			PosOut = list(range(LimInf,LimSup))
		elif LimSup >= gen_len:
			TurnSup = LimSup - gen_len
			PosOut = list(range(posMax,gen_len))+list(range(0,TurnSup)) + list(range(LimInf,posMax))
		elif LimInf < 0:
			TurnInf = gen_len + LimInf
			PosOut = list(range(0,posMax))+list(range(TurnInf,gen_len)) + list(range(posMax,LimSup))
		picMaxOK = []
		picOUT = []
		for peaks in picMaxRC:
			if peaks[1] not in PosOut:
				picMaxOK.append(peaks)
			else:
				picOUT.append(peaks)
		return picMaxOK, picOUT

	def addClosePic(picList, picClose, norm = 0):
		if norm:
			if picClose[0][0] >= 0.5:
				return picList, [picClose[0]]
		picListOK = picList[:]
		cov_add = 0
		for cov in picClose:
			cov_add += cov[0]
			picListOK[cov[1]] = 0.01
		picListOK[picClose[0][1]] = cov_add
		return picListOK, picClose

	def remove_pics(arr,n):
		arr=np.array(arr)
		pic_pos=arr.argsort()[-n:][::-1]
		arr2=np.delete(arr,pic_pos)
		return arr2

	def gamma(X):
		X = np.array(X, dtype=np.int64)
		v = remove_pics(X, 3)

		dist_max = float(max(v))
		if dist_max == 0:
			return np.array([1.00] * len(X))

		actual = np.bincount(v)
		fit_alpha, fit_loc, fit_beta = stats.gamma.fit(v)
		expected = stats.gamma.pdf(np.arange(0, dist_max + 1, 1), fit_alpha, loc=fit_loc, scale=fit_beta) * sum(actual)

		return stats.gamma.pdf(X, fit_alpha, loc=fit_loc, scale=fit_beta)


	def test_pics_decision_tree(whole_coverage, termini_coverage, termini_coverage_norm, termini_coverage_norm_close):
		L = len(whole_coverage[0])
		res = pd.DataFrame({"Position": np.array(range(L)) + 1, "termini_plus": termini_coverage[0],
                        "SPC_norm_plus": termini_coverage_norm[0], "SPC_norm_minus": termini_coverage_norm[1],
                        "SPC_norm_plus_close": termini_coverage_norm_close[0],
                        "SPC_norm_minus_close": termini_coverage_norm_close[1], "termini_minus": termini_coverage[1],
                        "cov_plus": whole_coverage[0], "cov_minus": whole_coverage[1]})

		res["cov"] = res["cov_plus"].values + res["cov_minus"].values

		res["R_plus"] = list(map(float, termini_coverage[0])) // np.mean(termini_coverage[0])
		res["R_minus"] = list(map(float, termini_coverage[1])) // np.mean(termini_coverage[1])

		regr = DecisionTreeRegressor(max_depth=3, min_samples_leaf=100)
		X = np.arange(L)
		X = X[:, np.newaxis]
		y = res["cov"].values
		regr.fit(X, y)

		y_1 = regr.predict(X)
		res["covnode"] = y_1
		covnodes = np.unique(y_1)
		thres = np.mean(whole_coverage[0]) / 2
		covnodes = [n for n in covnodes if n > thres]

		for node in covnodes:
			X = res[res["covnode"] == node]["termini_plus"].values
			res.loc[res["covnode"] == node, "pval_plus"] = gamma(X)
			X = res[res["covnode"] == node]["termini_minus"].values
			res.loc[res["covnode"] == node, "pval_minus"] = gamma(X)

		res.loc[res.pval_plus > 1, 'pval_plus'] = 1.00
		res.loc[res.pval_minus > 1, 'pval_minus'] = 1.00
		res = res.fillna(1.00)

		res['pval_plus_adj'] = multipletests(res["pval_plus"].values, alpha=0.01, method="bonferroni")[1]
		res['pval_minus_adj'] = multipletests(res["pval_minus"].values, alpha=0.01, method="bonferroni")[1]

		res = res.fillna(1.00)

		res_plus = pd.DataFrame(
			{"Position": res['Position'], "SPC_std": res['SPC_norm_plus'], "SPC": res['SPC_norm_plus_close'],
			"pval_gamma": res['pval_plus'], "pval_gamma_adj": res['pval_plus_adj']})
		res_minus = pd.DataFrame(
			{"Position": res['Position'], "SPC_std": res['SPC_norm_minus'], "SPC": res['SPC_norm_minus_close'],
			"pval_gamma": res['pval_minus'], "pval_gamma_adj": res['pval_minus_adj']})

		res_plus.sort_values("SPC", ascending=False, inplace=True)
		res_minus.sort_values("SPC", ascending=False, inplace=True)
	
		res_plus.reset_index(drop=True, inplace=True)
		res_minus.reset_index(drop=True, inplace=True)

		return res, res_plus, res_minus

	data = pd.read_csv("/home/emily/wf-NanoTerm/output/tau.csv")

	f = data[data['strand'] == "f"]
	r = data[data['strand'] == "r"]

	f_cov = f['cov'].to_numpy()
	r_cov = r['cov'].to_numpy()

	whole_coverage = np.array([f_cov, r_cov])

	f_spc = f['SPC'].to_numpy()
	r_spc = r['SPC'].to_numpy()

	termini_coverage = np.array([f_spc, r_spc])

	f_tau = f['tau'].to_numpy()
	r_tau = r['tau'].to_numpy()

	termini_coverage_norm = np.array([f_tau, r_tau])
	
	surrounding = 20
	gen_len = len(whole_coverage)

	picMaxPlus, picMaxMinus, TopFreqH = picMax(termini_coverage, 5)
	picMaxPlus_norm, picMaxMinus_norm, TopFreqH_norm = picMax(termini_coverage_norm, 5)

	picMaxPlus, picOUT_forw = RemoveClosePicMax(picMaxPlus, gen_len, surrounding)
	picMaxMinus, picOUT_rev = RemoveClosePicMax(picMaxMinus, gen_len, surrounding)
	picMaxPlus_norm, picOUT_norm_forw = RemoveClosePicMax(picMaxPlus_norm, gen_len, surrounding)
	picMaxMinus_norm, picOUT_norm_rev = RemoveClosePicMax(picMaxMinus_norm, gen_len, surrounding)

	termini_coverage_close = termini_coverage[:]
	termini_coverage_close[0], picOUT_forw = addClosePic(termini_coverage[0], picOUT_forw)
	termini_coverage_close[1], picOUT_rev = addClosePic(termini_coverage[1], picOUT_rev)

	termini_coverage_norm_close = termini_coverage_norm[:]
	termini_coverage_norm_close[0], picOUT_norm_forw = addClosePic(termini_coverage_norm[0], picOUT_norm_forw, 1)
	termini_coverage_norm_close[1], picOUT_norm_rev = addClosePic(termini_coverage_norm[1], picOUT_norm_rev, 1)

	phage_norm, phage_plus_norm, phage_minus_norm = test_pics_decision_tree(whole_coverage, termini_coverage, termini_coverage_norm, termini_coverage_norm_close)

	np.savetxt("phage_norm.csv", phage_norm, delimiter=",")
	"""
}

process extMapping {
	input:
		path input
		path extref
		val seqplat

	output:
		path 'ext_aln.sam'
		
	script:
		if( seqplat == 'nanopore' )
			"""
			minimap2 -ax map-ont $extref $params.input > ext_aln.sam
			"""
		
		else if( seqplat == 'illumina' )
			"""
			minimap2 -asr $extref $params.input > ext_aln.sam
			"""

		else 
			 error "Invalid sequencing platform: $params.seqplat"
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
	samtools depth -a -d 0 -o ext_aln_f_cov.txt $ext_aln_f_sorted
	samtools depth -a -d 0 -o ext_aln_r_cov.txt $ext_aln_r_sorted
	"""
}

process extTau {
	publishDir "${params.outdir}", mode: 'copy', overwrite: true
	
	input:
		path ext_aln_f
		path ext_aln_r
		path ext_aln_f_cov
		path ext_aln_r_cov
		path extFasta
		
	output:
		path 'trim_tau.csv'
		
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
                    colClasses=c("character","numeric","numeric"))
	f_cov <- setNames(f_cov, c("chr","pos","cov"))
	f_cov['strand'] <- "f"

	r_cov <- read.table("$ext_aln_r_cov", quote="\\"", comment.char="",
                    colClasses=c("character","numeric","numeric"))
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
	
	trim_tau <- subset(ext_tau, ext_tau['pos'] > 500)
	trim_tau['pos'] <- trim_tau['pos'] - 500
	max_pos <- len - 1000
	trim_tau <- subset(trim_tau, trim_tau['pos'] <= max_pos)
	
	write.csv(trim_tau, "trim_tau.csv")
	"""	
}

workflow {
	len_ch = rawSeq(ref_ch)
	ext_ch = extendRef(len_ch)
	aln_ch = mapping(seq_ch, ref_ch, plat_ch)
	sep_ch = strandSep(aln_ch)
	bed_ch = bed(sep_ch)
	cov_ch = cov(sep_ch)
	tau_ch = tau(bed_ch, cov_ch, len_ch)
	stats(tau_ch)
	extlen_ch = extRawSeq(ext_ch)
	extaln_ch = extMapping(seq_ch, ext_ch, plat_ch)
	extsep_ch = extStrandSep(extaln_ch)
	extbed_ch = extBed(extsep_ch)
	extcov_ch = extCov(extsep_ch)
	extTau(extbed_ch, extcov_ch, extlen_ch)
}