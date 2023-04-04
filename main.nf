#!/usr/bin/env nextflow

ref_ch = Channel.of(params.fasta)
seq_ch = Channel.of(params.input)
plat_ch = Channel.of(params.seqplat)
name_ch = Channel.of(params.name)
outdir_ch = Channel.of(params.outdir)

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
	publishDir "${params.outdir}", mode: 'copy', overwrite: true
	
	input:
		path input
		path fasta
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
		body_add_par(value = paste("Input sequences: ", "$params.input"), style = "Normal") %>%
		body_add_par(value = paste("Reference genome: ", "$params.fasta"), style = "Normal") %>%
		body_add_par(value = paste("Reference genome length: ", len), style = "Normal") %>%
		body_add_par("Alignment details", style = "heading 2") %>%
		body_add_par(value = paste("Sequencing platform:", "$params.seqplat"), style = "Normal") %>%
		body_add_par(value = paste("Number of sequence reads: ", totalReads), style = "Normal") %>%
		body_add_par(value = paste("Number of reads aligned: ", mappedReads), style = "Normal") %>%
		body_add_par(value = paste("Number of reads not aligned: ", unmappedReads), style = "Normal") %>%
		body_add_par(value = paste("Average read length: ", aveReadLen), style = "Normal") %>%
		body_add_par(value = paste("Maximum read length: ", maxReadLen), style = "Normal") %>%
		body_add_par(value = paste("Average read depth: ", meanDepth), style = "Normal") %>%
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
	publishDir "${params.outdir}", mode: 'copy', overwrite: true
	
	input:
		path report
		path outdir
		
	output:
		path 'report.pdf', optional: true
	
	script:
	"""
	libreoffice --headless --convert-to pdf  $report --outdir $params.outdir
	"""	
}

workflow {
	len_ch = rawSeq(ref_ch)
	aln_ch = mapping(seq_ch, ref_ch, plat_ch)
	alnstats_ch = alignStats(aln_ch)
	sep_ch = strandSep(aln_ch)
	bed_ch = bed(sep_ch)
	cov_ch = cov(sep_ch)
	tau_ch = tau(bed_ch, cov_ch, len_ch)
	stats_ch = stats(tau_ch)
	report_ch = report(stats_ch, ref_ch, seq_ch, name_ch, plat_ch, alnstats_ch)
	doc2pdf(report_ch, outdir_ch)
}