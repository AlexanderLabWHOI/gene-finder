configfile: "config.yaml"

import io
import os
import glob
import pandas as pd
import numpy as np
import pathlib
from snakemake.exceptions import print_exception, WorkflowError

#----SET VARIABLES----#
PROTEIN_DIR = config['protein_dir']
ALIGNMENT_DIR = config['alignment_dir']
MAG_DIR = config['mag_dir']
OUTPUT_DIR = config['output_dir']
GENEFILE = config['gene_list']
genelist = []
with open(GENEFILE, 'r') as f:
    for line in f:
        genelist.append(line.strip())

SAMPLES = [os.path.basename(f).replace(".proteins.faa", "") for f in glob.glob(PROTEIN_DIR + "/*.proteins.faa")]

#----RULES----#
localrules: parse_hmmsearch

rule all:
    input:
        hmmbuild =  expand('{base}/{gene}/gene_alignment_profile.hmm', base = ALIGNMENT_DIR, gene=genelist), 
        hmmsearch = expand('{base}/hmm_results/{gene}/{sample}.hmmout', base = OUTPUT_DIR, gene = genelist, sample = SAMPLES ), 
        hmmtable = expand('{base}/hmm_results/{gene}/{sample}.eval.tab', base = OUTPUT_DIR, gene = genelist, sample = SAMPLES ), 
        hmm_allhits = expand('{base}/hmm_results/{gene}/{sample}.hits.faa', base = OUTPUT_DIR, gene = genelist, sample = SAMPLES ),
        hmm_mags = expand("{base}/mag_results/{gene}.maghits.contigs.csv", base = OUTPUT_DIR, gene = genelist) 
 

rule hmmbuild:
    input: alignment = ALIGNMENT_DIR + "/{gene}/gene_alignment.aln"
    output: hmm = ALIGNMENT_DIR + "/{gene}/gene_alignment_profile.hmm"
    conda: 
        "envs/hmmer.yaml"
    shell:
        """
        hmmbuild {output.hmm} {input.alignment} 
        """

rule hmmsearch:
    input: 
        proteins = PROTEIN_DIR + "/{sample}.proteins.faa", 
        hmm = ALIGNMENT_DIR + "/{gene}/gene_alignment_profile.hmm"
    output: 
        hmmout = OUTPUT_DIR + "/hmm_results/{gene}/{sample}.hmmout", 
    params:
        cpu = "--cpu 2"
    conda:
        "envs/hmmer.yaml"
    shell:
        """
        hmmsearch {params.cpu} {input.hmm} {input.proteins} > {output.hmmout}  
        """

rule parse_hmmsearch:
    input:
        OUTPUT_DIR + "/hmm_results/{gene}/{sample}.hmmout"
    output:
        OUTPUT_DIR + "/hmm_results/{gene}/{sample}.eval.tab"
    conda:
        "envs/biopython.yaml"
    script:
        "scripts/bioconda-parser.py"

rule get_contig_hits:
    input:
        eval_tab = OUTPUT_DIR + "/hmm_results/{gene}/{sample}.eval.tab", 
        proteins = PROTEIN_DIR + "/{sample}.proteins.faa",
    output:
        OUTPUT_DIR + "/hmm_results/{gene}/{sample}.hits.faa"
    conda:
        "envs/seqtk.yaml"
    shell: 
        """
        cut -f1 {input.eval_tab} | seqtk subseq {input.proteins} - > {output}
        """

rule get_MAG_hits:
    input:
        genelist = OUTPUT_DIR + "/curated-gene-lists/{gene}.threshold.csv", 
        magIDList = MAG_DIR + "/mag_fasta_headers.txt"
    output:
        maglist = OUTPUT_DIR + "/mag_results/{gene}.maghits.contigs.csv"
    conda: "envs/biopython.yaml" 
    shell: 
        """
        python scripts/mag-finder.py {input.genelist} {input.magIDList} {output.maglist}        
        """
