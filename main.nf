#! /usr/bin/env nextflow

blastdb="MyBlastDatabase"
params.query="file.fasta"

println "I will BLAST $params.query against $blastdb"

