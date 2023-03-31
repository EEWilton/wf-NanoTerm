#!/usr/bin/env nextflow

ref_ch = Channel.of(params.fasta)
seq_ch = Channel.of(params.input)
plat_ch = Channel.of(params.seqplat)
name_ch = Channel.of(params.name)

println "\nReference: $params.fasta\nSequence reads: $params.input\n"

// This process extracts the raw sequence from the fasta file, removing the header
// The output is a text file with just the nucleotide sequence of the reference genome
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

// This process extends the raw referemce sequence by 500 nt on each end, as if the genome was a concatamer
// The output is a new fasta file with the extended reference sequence
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

// This process extracts the extended reference sequence into a text file, removing the header
// The output is a text file with just the nucleotide sequence of the extended reference genome
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

// This process maps the input fastq.gz file against the reference genome using minimap2
// This process can be adjusted to suit the sequencing platform used
// The output is a .sam file of the alignment
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

process alignStats {
	input:
		path aln
		
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
		path 'stats.csv'
		path 'f_sig.csv'
		path 'r_sig.csv'
	
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

	def peakMax(coverage, nbr_peaks):
		peakMaxPlus = heapq.nlargest(nbr_peaks, zip(coverage[0], itertools.count()))
		peakMaxMinus = heapq.nlargest(nbr_peaks, zip(coverage[1], itertools.count()))
		TopFreqH = max(max(np.array(list(zip(*peakMaxPlus))[0])), max(np.array(list(zip(*peakMaxMinus))[0])))
		return peakMaxPlus, peakMaxMinus, TopFreqH

	def removeClosePeakMax(peakMax, gen_len, nbr_base):
		if nbr_base == 0:
			return peakMax[1:], [peakMax[0]]
		peakMaxRC = peakMax[:]
		posMax = peakMaxRC[0][1]
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
		peakMaxOK = []
		peakOUT = []
		for peaks in peakMaxRC:
			if peaks[1] not in PosOut:
				peakMaxOK.append(peaks)
			else:
				peakOUT.append(peaks)
		return peakMaxOK, peakOUT

	def addClosePeak(peakList, peakClose, norm = 0):
		if norm:
			if peakClose[0][0] >= 0.5:
				return peakList, [peakClose[0]]
		peakListOK = peakList[:]
		cov_add = 0
		for cov in peakClose:
			cov_add += cov[0]
			peakListOK[cov[1]] = 0.01
		peakListOK[peakClose[0][1]] = cov_add
		return peakListOK, peakClose

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

	def peaksDecisionTree(read_depth, starting_pos_depth, tau, tau_close):
		L = len(read_depth[0])
		res = pd.DataFrame({"Position": np.array(range(L)) + 1, "SPC_plus": starting_pos_depth[0],
                        "tau_plus": tau[0], "tau_minus": tau[1],
                        "tau_plus_close": tau_close[0],
                        "tau_minus_close": tau_close[1], "SPC_minus": starting_pos_depth[1],
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
			{"Position": res['Position'], "tau": res['tau_plus'], "tau_close": res['tau_plus_close'],
			"pval_gamma": res['pval_plus'], "pval_gamma_adj": res['pval_plus_adj']})
		res_minus = pd.DataFrame(
			{"Position": res['Position'], "tau": res['tau_minus'], "tau_close": res['tau_minus_close'],
			"pval_gamma": res['pval_minus'], "pval_gamma_adj": res['pval_minus_adj']})

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

	peakMaxPlus, peakMaxMinus, TopFreqH = peakMax(starting_pos_depth, 5)
	peakMaxPlus_norm, peakMaxMinus_norm, TopFreqH_norm = peakMax(tau, 5)

	peakMaxPlus, peakOUT_forw = removeClosePeakMax(peakMaxPlus, gen_len, surrounding)
	peakMaxMinus, peakOUT_rev = removeClosePeakMax(peakMaxMinus, gen_len, surrounding)
	peakMaxPlus_norm, peakOUT_norm_forw = removeClosePeakMax(peakMaxPlus_norm, gen_len, surrounding)
	peakMaxMinus_norm, peakOUT_norm_rev = removeClosePeakMax(peakMaxMinus_norm, gen_len, surrounding)

	starting_pos_depth_close = starting_pos_depth[:]
	starting_pos_depth_close[0], peakOUT_forw = addClosePeak(starting_pos_depth[0], peakOUT_forw)
	starting_pos_depth_close[1], peakOUT_rev = addClosePeak(starting_pos_depth[1], peakOUT_rev)

	tau_close = tau[:]
	tau_close[0], peakOUT_norm_forw = addClosePeak(tau[0], peakOUT_norm_forw, 1)
	tau_close[1], peakOUT_norm_rev = addClosePeak(tau[1], peakOUT_norm_rev, 1)

	phage_norm, phage_plus_norm, phage_minus_norm = peaksDecisionTree(read_depth, starting_pos_depth, tau, tau_close)
	
	plus_significant = selectSignificant(phage_plus_norm, 1.0 / gen_len, 1.0)
	minus_significant = selectSignificant(phage_minus_norm, 1.0 / gen_len, 1.0)

	plus_significant.to_csv("f_sig.csv")
	minus_significant.to_csv("r_sig.csv")
	phage_norm.to_csv("stats.csv")
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

process extStats {
	publishDir "${params.outdir}", mode: 'copy', overwrite: true

	input:
		path trim_tau
		
	output:
		path 'ext_stats.csv'
		path 'ext_f_sig.csv'
		path 'ext_r_sig.csv'
	
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

	def peakMax(coverage, nbr_peaks):
		peakMaxPlus = heapq.nlargest(nbr_peaks, zip(coverage[0], itertools.count()))
		peakMaxMinus = heapq.nlargest(nbr_peaks, zip(coverage[1], itertools.count()))
		TopFreqH = max(max(np.array(list(zip(*peakMaxPlus))[0])), max(np.array(list(zip(*peakMaxMinus))[0])))
		return peakMaxPlus, peakMaxMinus, TopFreqH

	def removeClosePeakMax(peakMax, gen_len, nbr_base):
		if nbr_base == 0:
			return peakMax[1:], [peakMax[0]]
		peakMaxRC = peakMax[:]
		posMax = peakMaxRC[0][1]
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
		peakMaxOK = []
		peakOUT = []
		for peaks in peakMaxRC:
			if peaks[1] not in PosOut:
				peakMaxOK.append(peaks)
			else:
				peakOUT.append(peaks)
		return peakMaxOK, peakOUT

	def addClosePeak(peakList, peakClose, norm = 0):
		if norm:
			if peakClose[0][0] >= 0.5:
				return peakList, [peakClose[0]]
		peakListOK = peakList[:]
		cov_add = 0
		for cov in peakClose:
			cov_add += cov[0]
			peakListOK[cov[1]] = 0.01
		peakListOK[peakClose[0][1]] = cov_add
		return peakListOK, peakClose

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


	def peaksDecisionTree(read_depth, starting_pos_depth, tau, tau_close):
		L = len(read_depth[0])
		res = pd.DataFrame({"Position": np.array(range(L)) + 1, "SPC_plus": starting_pos_depth[0],
                        "tau_plus": tau[0], "tau_minus": tau[1],
                        "tau_plus_close": tau_close[0],
                        "tau_minus_close": tau_close[1], "SPC_minus": starting_pos_depth[1],
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
			{"Position": res['Position'], "tau": res['tau_plus'], "tau_close": res['tau_plus_close'],
			"pval_gamma": res['pval_plus'], "pval_gamma_adj": res['pval_plus_adj']})
		res_minus = pd.DataFrame(
			{"Position": res['Position'], "tau": res['tau_minus'], "tau_close": res['tau_minus_close'],
			"pval_gamma": res['pval_minus'], "pval_gamma_adj": res['pval_minus_adj']})

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

	data = pd.read_csv("$trim_tau")

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
	gen_len = len(read_depth)

	peakMaxPlus, peakMaxMinus, TopFreqH = peakMax(starting_pos_depth, 5)
	peakMaxPlus_norm, peakMaxMinus_norm, TopFreqH_norm = peakMax(tau, 5)

	peakMaxPlus, peakOUT_forw = removeClosePeakMax(peakMaxPlus, gen_len, surrounding)
	peakMaxMinus, peakOUT_rev = removeClosePeakMax(peakMaxMinus, gen_len, surrounding)
	peakMaxPlus_norm, peakOUT_norm_forw = removeClosePeakMax(peakMaxPlus_norm, gen_len, surrounding)
	peakMaxMinus_norm, peakOUT_norm_rev = removeClosePeakMax(peakMaxMinus_norm, gen_len, surrounding)

	starting_pos_depth_close = starting_pos_depth[:]
	starting_pos_depth_close[0], peakOUT_forw = addClosePeak(starting_pos_depth[0], peakOUT_forw)
	starting_pos_depth_close[1], peakOUT_rev = addClosePeak(starting_pos_depth[1], peakOUT_rev)

	tau_close = tau[:]
	tau_close[0], peakOUT_norm_forw = addClosePeak(tau[0], peakOUT_norm_forw, 1)
	tau_close[1], peakOUT_norm_rev = addClosePeak(tau[1], peakOUT_norm_rev, 1)

	phage_norm, phage_plus_norm, phage_minus_norm = peaksDecisionTree(read_depth, starting_pos_depth, tau, tau_close)
	
	plus_significant = selectSignificant(phage_plus_norm, 1.0 / gen_len, 1.0)
	minus_significant = selectSignificant(phage_minus_norm, 1.0 / gen_len, 1.0)

	plus_significant.to_csv("ext_f_sig.csv")
	minus_significant.to_csv("ext_r_sig.csv")

	phage_norm.to_csv("ext_stats.csv")
	"""
}

process report {
	publishDir "${params.outdir}", mode: 'copy', overwrite: true
	
	input:
		path input
		path fasta
		path stats
		path f_sig
		path r_sig
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
		arrange(desc(tau_plus)) %>%
		filter(pval_plus_adj <= 0.5) %>%
		filter(cov_plus > (0.05 * max_plus)) 
	
	top_tau_minus <- data %>%
		arrange(desc(tau_minus)) %>%
		filter(pval_minus_adj <= 0.5) %>%
		filter(cov_minus > (0.05 * max_minus)) 

	f_term <- top_tau_plus[1,2]
	r_term <- top_tau_minus[1,2]

	f_term_tau <- top_tau_plus[1,4]
	r_term_tau <- top_tau_minus[1,5]
	
	f_term_p <- top_tau_plus[1,17]
	r_term_p <- top_tau_minus[1,18]

	strand <- c("Forward","Reverse")
	terms <- c(f_term, r_term)
	taus <- c(f_term_tau, r_term_tau)
	pvals <- c(f_term_p, r_term_p)
	table <- data.frame(strand, terms, taus, pvals)

	term_dist <- r_term - f_term
	
	colours <- c("Forward" = "springgreen4", "Reverse" = "purple")

	depth <- ggplot(data = data, aes(x=Position, y=rollmean(cov, window, na.pad = TRUE, align = "right"))) +
		theme_calc() + 
		geom_vline(xintercept=f_term, linetype="dashed", colour="springgreen3", linewidth=1.1) + 
		geom_vline(xintercept=r_term, linetype="dashed", colour="violet", linewidth=1.1) +
		geom_line(data=data, aes(x=Position, y=rollmean(cov_plus, window, na.pad = TRUE, align = "right"), colour="Forward")) +
		geom_line(data=data, aes(x=Position, y=rollmean(cov_minus, window, na.pad = TRUE, align = "right"), colour="Reverse")) +
		geom_line() +
		labs(x = "Reference genome position",
			y = "Read depth",
			colour = "Legend") +
		scale_color_manual(values = colours) +
		scale_x_continuous(labels = comma) +
		guides(colour = guide_legend(override.aes = list(linewidth = 3)))

	tau <- ggplot(data=data) +
		theme_calc() + 
		geom_point(data=data, aes(x=Position, y=tau_plus, colour="Forward")) +
		geom_point(data=data, aes(x=Position, y=tau_minus, colour="Reverse")) +
		geom_label_repel(data=subset(data, Position == f_term), 
                   aes(x=Position, y=tau_plus,label=Position), colour="springgreen4",
                   show.legend = FALSE) + 
		geom_label_repel(data=subset(data,  Position == r_term), 
                   aes(x=Position, y=tau_minus,label=Position), colour="purple",
                   show.legend = FALSE) +
		labs(x = "Reference genome position",
			y = "tau",
			colour = "Legend") +
		scale_color_manual(values = colours) +
		scale_x_continuous(labels = comma) +
		scale_y_continuous(limits=c(0,1))
	
	
	if (term_dist > 20){
		class <- "DTR"
	} else if (term_dist < 0){
		class <- "COS 3′"
	} else {
		class <- "COS 5′"
	}

	report <- read_docx() %>%
		body_add_par(value = paste("NanoTerm Report: ", name), style = "heading 1") %>%
		body_add_par("", style = "Normal") %>%
		body_add_par("Run details", style = "heading 2") %>%
		body_add_par(value = paste("Generated on: ", date), style = "Normal") %>%
		body_add_par(value = paste("Input sequences: ", "$params.input"), style = "Normal") %>%
		body_add_par(value = paste("Reference genome: ", "$params.fasta"), style = "Normal") %>%
		body_add_par("Alignment details", style = "heading 2") %>%
		body_add_par(value = paste("Sequencing platform:", "$params.seqplat"), style = "Normal") %>%
		body_add_par(value = paste("Number of sequence reads: ", totalReads), style = "Normal") %>%
		body_add_par(value = paste("Number of reads aligned: ", mappedReads), style = "Normal") %>%
		body_add_par(value = paste("Number of reads not aligned: ", unmappedReads), style = "Normal") %>%
		body_add_par(value = paste("Average read length: ", aveReadLen), style = "Normal") %>%
		body_add_par(value = paste("Maximum read length: ", maxReadLen), style = "Normal") %>%
		body_add_par("", style = "Normal") %>%
		body_add_par("Phage prediction", style = "heading 2") %>%
		body_add_table(table, style = "Normal", first_column = TRUE) %>%
		body_add_par("", style = "Normal") %>%
		body_add_par(value = paste("Phage class: ", class), style = "Normal")

	if (class == "COS 5′"){
		report <- body_add_par(report, value = paste("Cohesive sequence: ", class), style = "Normal")
	} else if (class == "COS 3′"){
		report <- body_add_par(report, value = paste("Cohesive sequence: ", class), style = "Normal")
	} else if (class == "DTR"){
		report <- body_add_par(report, value = paste("DTR length: ", term_dist), style = "Normal")
	} else {
	}

	report <- body_add_par(report, "Figures", style = "heading 2") %>%
		body_add_par("", style = "Normal") %>%
		body_add_gg(value = tau, style = "centered") %>%
		body_add_par(value = "Figure 1. The tau value calculated for each genome position.") %>%
		body_add_par("", style = "Normal") %>%
		body_add_gg(value = depth, style = "centered") %>%
		body_add_par(value = "Figure 2. The total read depth of the sequencing run, graphed as a rolling average with a window size equal to 1% of the reference genome length.  Black is the sum of forward and reverse read depth.")
 
	print(report, target = "./report.docx")
	"""
}

workflow {
	len_ch = rawSeq(ref_ch)
	ext_ch = extendRef(len_ch)
	aln_ch = mapping(seq_ch, ref_ch, plat_ch)
	alnstats_ch = alignStats(aln_ch)
	sep_ch = strandSep(aln_ch)
	bed_ch = bed(sep_ch)
	cov_ch = cov(sep_ch)
	tau_ch = tau(bed_ch, cov_ch, len_ch)
	stats_ch = stats(tau_ch)
	extlen_ch = extRawSeq(ext_ch)
	extaln_ch = extMapping(seq_ch, ext_ch, plat_ch)
	extsep_ch = extStrandSep(extaln_ch)
	extbed_ch = extBed(extsep_ch)
	extcov_ch = extCov(extsep_ch)
	extTau_ch = extTau(extbed_ch, extcov_ch, extlen_ch)
	extStats_ch = extStats(extTau_ch)
	report(stats_ch, ref_ch, seq_ch, name_ch, plat_ch, alnstats_ch)
}
