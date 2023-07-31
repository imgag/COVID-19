import pathlib
import re
from datetime import datetime
import os
import glob

# Load Configfile in Repo and WorkingDir
configfile: srcdir('config.yaml')
if os.path.isfile('config.yaml'):
  print("Found project config file")
  configfile: 'config.yaml'

def get_samples():
    if config["samples"] is not None:
      return [item.strip().replace("Sample_", "") for item in config["samples"]]
    else:
      return [x.replace("Sample_", "") for x in glob.glob("Sample_*") if os.path.isdir(x)]

def get_fastqs_by_R1(wildcards):
    fastqs = glob.glob("Sample_%s/%s*L00*R1*.fastq.gz"%(wildcards.s,wildcards.s))
    return sorted(fastqs,key=lambda x: x)
    
def get_fastqs_by_R2(wildcards):
    fastqs = glob.glob("Sample_%s/%s*L00*R2*.fastq.gz"%(wildcards.s,wildcards.s))
    return sorted(fastqs,key=lambda x: x)

def get_fastqs_by_index(wildcards):
    fastqs = glob.glob("Sample_%s/%s*index*.fastq.gz"%(wildcards.s,wildcards.s))
    return sorted(fastqs,key=lambda x: x)

all_sampleids = get_samples()
threads_max = config['threads']

rule all:
  input:
    config['krakendb'], # krakendb 
    config['reference'],    # reference genome
    expand('Sample_{s}/{s}.viral.bam', s = all_sampleids), # mapping
    expand('Sample_{s}/coverage/{s}.mosdepth.summary.txt', s = all_sampleids), # coverage depth qc
    expand('Sample_{s}/ivar/{s}_ivar.vcf.gz', s = all_sampleids), # variant calling
    expand('Sample_{s}/consensus/{s}_consensus_ivar.fa', s = all_sampleids),  # consensus ivar
    #expand('Sample_{s}/consensus/{s}_umivar.fasta', s = all_sampleids),  # consensus umivar
    expand('Sample_{s}/quast_results', s = all_sampleids), # assembly 
    expand('Sample_{s}/ivar/{s}.variants.txt', s = all_sampleids), # variants.txt
    #expand('Sample_{s}/umivar2/{s}.viral_hq.vcf', s = all_sampleids), #umiVar2
    expand('Sample_{s}/umivar2/{s}_bamclipoverlap_sorted_hq.vcf', s = all_sampleids), #umiVar2
    expand('Sample_{s}/lofreq/{s}_lofreq.tsv', s = all_sampleids), # Lofreq
    expand('Sample_{s}/varscan/{s}_varscan.tsv', s = all_sampleids), # varscan
    expand('Sample_{s}/dedup/ivar/{s}.variants.txt', s = all_sampleids), # variants.txt
    expand('Sample_{s}/dedup/lofreq/{s}_lofreq.tsv', s = all_sampleids), # Lofreq
    expand('Sample_{s}/dedup/varscan/{s}_varscan.tsv', s = all_sampleids), # varscan
    expand('Sample_{s}/{s}_R1_trimmed.fastq.gz', s = all_sampleids)#want to check output 


##why here?
# INDEX REFERENCE genome
rule index_reference:
  input: 
     ref= config['reference'] 
  output:
    amb = '%s.amb'%config['reference'],
    ann = '%s.ann'%config['reference'],
    bwt = '%s.bwt'%config['reference'],
    pac = '%s.pac'%config['reference'],
    sa = '%s.sa'%config['reference']
  conda: 'envs/env_bwa.yaml'
  shell:
    """
    bwa index {input.ref}
    """
		
#_______ READ PREPROCESSING__________________________________________________________# 

# Add UMI Barcodes to FASTQ header
rule add_barcode:
    input:
        in1=get_fastqs_by_R1,
        in2=get_fastqs_by_R2,
        barcode=get_fastqs_by_index
    output:
        #out1=temp('Sample_{s}/{s}_R1_barcode_added.fastq.gz'),
        #out2=temp('Sample_{s}/{s}_R2_barcode_added.fastq.gz')
        out1='Sample_{s}/{s}_R1_barcode_added.fastq.gz',
        out2='Sample_{s}/{s}_R2_barcode_added.fastq.gz'
    conda: 'envs/env_ngs_bits.yaml'
    shell:
        """
        FastqAddBarcode -in1 {input.in1} -in2 {input.in2} -in_barcode {input.barcode} -out1 {output.out1} -out2 {output.out2}
        """

# Read trimming and ReadQC
rule trimming:
    input:
        in1=rules.add_barcode.output.out1,
        in2=rules.add_barcode.output.out2
    output:
        #out1=temp('Sample_{s}/{s}_R1_trimmed.fastq.gz'),#want to keep temp files to test newer version of kraken
        #out2=temp('Sample_{s}/{s}_R2_trimmed.fastq.gz'),
        out1='Sample_{s}/{s}_R1_trimmed.fastq.gz',
        out2='Sample_{s}/{s}_R2_trimmed.fastq.gz',
        qc='Sample_{s}/{s}_stats_fastq.qcML'
    threads: 
        threads_max
    conda: 'envs/env_ngs_bits.yaml'
    params:
        a1=config['adapter1'],
	a2=config['adapter2']
    shell:
        """
        SeqPurge -in1 {input.in1} -in2 {input.in2} -out1 {output.out1} -out2 {output.out2} -a1 {params.a1} -a2 {params.a2} -qc {output.qc} -threads {threads}
        """

# Removal of human host reads with Kraken 
rule kraken:
    input:
        in1=rules.trimming.output.out1,
        in2=rules.trimming.output.out2
    output:
        out1='Sample_{s}/{s}_viral_1.fastq.gz',
        out2='Sample_{s}/{s}_viral_2.fastq.gz',
    threads: 
        threads_max
    conda: 'envs/env_kraken.yaml'
    params:
        db=config['krakendb']
    shell:
       """
       kraken2  \
        --db {params.db} \
        -threads {threads} \
        --unclassified-out Sample_{wildcards.s}/{wildcards.s}_viral#.fastq  \
        --report Sample_{wildcards.s}/{wildcards.s}.kraken2.report.txt  \
        --report-zero-counts  \
        --paired --gzip-compressed \
        {input.in1} {input.in2} > /dev/null
       gzip Sample_{wildcards.s}/{wildcards.s}_viral_1.fastq
       gzip Sample_{wildcards.s}/{wildcards.s}_viral_2.fastq
       """


#_______ MAPPING ___________________________________________________________#

rule mapping:
  input: 
     r1=rules.kraken.output.out1,
     r2=rules.kraken.output.out2,
     index_files = rules.index_reference.output
  output:
    bam = 'Sample_{s}/{s}.viral.bam',
    idx = 'Sample_{s}/{s}.viral.bam.bai'
  threads: threads_max
  params: 
    ref = config['reference']
  conda: 'envs/env_bwa.yaml'
  shell:
    """
    bwa mem -t {threads} {params.ref} {input.r1} {input.r2} -M | \
      samtools view -bS - | \
      samtools sort -@4 -m 1G -o {output.bam}
    samtools index {output.bam}
    """
		
rule mapping_qc:
  input:
    bam = rules.mapping.output.bam
  output:
    qcml = 'Sample_{s}/{s}_stats_map.qcML'
  params:
    target = config["target_twist"]
  conda: 'envs/env_ngs_bits.yaml'
  shell:
    """
    MappingQC -out {output.qcml} -roi {params.target} -no_cont -in {input.bam}
    """

rule cov_depth_qc:
  input:
    bam = rules.mapping.output.bam
  output:
    qcml = 'Sample_{s}/coverage/{s}.mosdepth.summary.txt'
  params:
    target = config["target_twist"]
  conda: 'envs/env_mosdepth.yaml'
  threads: threads_max
  shell:
    """
    mosdepth --by {params.target} -t {threads} Sample_{wildcards.s}/coverage/{wildcards.s} {input}
    """
#duplicate samblaster
rule mapping_duplicate_samblaster:
  input: 
    bam=rules.mapping.output.bam
  output:
    bam = 'Sample_{s}/dedup/{s}.viral.samblaster.bam',
    idx = 'Sample_{s}/dedup/{s}.viral.samblaster.bam.bai'
  threads: threads_max
  params: 
    ref = config['reference']
  conda: 'envs/env_bwa.yaml'
  shell:
    """

    samtools sort -@4 -m 1G -n {input.bam}|samtools view -h | samblaster -M | samtools sort -o {output.bam}
    samtools index {output.bam}
    """

#_______ DEDUPLICATE UMI________________________________________________________#

rule barcode_correction:
  input:
    bam =  rules.mapping.output.bam
  output: 
    bam = 'Sample_{s}/dedup/{s}.viral.corrected.bam'
  conda:
     "envs/env_umivar.yaml"
  log:
    "logs/{s}_barcode_correction.log"    
  params:
    barcode_correction = '%s/barcode_correction.py'%config['umiVar']
  shell:
    """
    {params.barcode_correction} --infile {input.bam} --outfile {output.bam} > {log}
    """

rule sort_by_position_after_correction:
    input: 
        rules.barcode_correction.output.bam
    output:
        corrected_sorted_bam = 'Sample_{s}/dedup/{s}_corrected_sorted.bam'
    conda: 'envs/env_samtools.yaml'
    threads: threads_max
    shell:"samtools sort -@ {threads} -m 1G -o {output.corrected_sorted_bam} {input}"	

rule index_after_correction:
    input: 
        corrected_sorted_bam = rules.sort_by_position_after_correction.output.corrected_sorted_bam
    output:
        'Sample_{s}/dedup/{s}_corrected_sorted.bam.bai'
    conda: 'envs/env_samtools.yaml'
    shell:
        "samtools index {input.corrected_sorted_bam}"	

rule bamclipoverlap:
    input:
       corrected_sorted_bam = rules.sort_by_position_after_correction.output.corrected_sorted_bam,
       corrected_sorted_idx = rules.index_after_correction.output
    output:
       overlapped_bam = 'Sample_{s}/dedup/{s}_bamclipoverlap.bam'
    conda: 'envs/env_ngs_bits.yaml'
    shell: 
       """
	   BamClipOverlap -in {input.corrected_sorted_bam} -out {output.overlapped_bam} -overlap_mismatch_basen
       """

rule sort_by_position_overlap:
    input: 
        rules.bamclipoverlap.output
    output:
        bamclipoverlap_sorted_bam = 'Sample_{s}/dedup/{s}_bamclipoverlap_sorted.bam'
    conda: 'envs/env_samtools.yaml'
    threads: threads_max
    shell:"samtools sort -@ {threads} -m 1G -o {output.bamclipoverlap_sorted_bam} {input}"	

rule index_overlap:
    input: 
        bamclipoverlap_sorted_bam = rules.sort_by_position_overlap.output.bamclipoverlap_sorted_bam
    output:
        idx = 'Sample_{s}/dedup/{s}_bamclipoverlap_sorted.bam.bai'
    conda: 'envs/env_samtools.yaml'
    shell:
        "samtools index {input.bamclipoverlap_sorted_bam}"	

#_______ VARIANT CALLING ________________________________________________________#

# Umivar2
rule umivar2:
  input: 
    bam = rules.sort_by_position_overlap.output.bamclipoverlap_sorted_bam#rules.mapping.output.bam#'Sample_{s}/{s}.viral.bam'
  output: 'Sample_{s}/umivar2/{s}_bamclipoverlap_sorted_hq.vcf'#'Sample_{s}/umivar2/{s}.viral_hq.vcf'
  conda: 'envs/env_umivar.yaml'
  log:
      "logs/{s}_umivar2.log"
  params:
    target = config["target_twist"],
    ref = config['reference'],
    umiVar= '%s/umiVar.py'%config['umiVar'],
    ac = "-ac " + config['umivar']['ac'] if 'ac' in config['umivar'] else "",
    af = "-af " + config['umivar']['af'] if 'af' in config['umivar'] else "",
    ns = "-ns " + config['umivar']['ns'] if 'ns' in config['umivar'] else "",
    sb = "-sb " + config['umivar']['sb'] if 'sb' in config['umivar'] else "",
    kt = "-kt ",# if 'kt' in config['umivar'] else "",
    mq = "-mq 30", #default
    bq = "-bq 20"  #default
  shell:
    """
    {params.umiVar} \
        -tbam {input.bam} \
        -b {params.target} \
        -r {params.ref} \
        -o Sample_{wildcards.s}/umivar2 \
        {params.ac} {params.af} {params.ns} {params.sb} {params.kt}  {params.mq}  {params.bq} 
    """
#umivar:  ac: "4"  ns: "-1" sb: "0" mq "30" bq "20"

# LoFreq
rule lofreq_call:
    input:
        bam = rules.sort_by_position_overlap.output.bamclipoverlap_sorted_bam
    output:
        "Sample_{s}/lofreq/{s}_lofreq.tsv"
    params:
        ref = config['reference']
    conda: 'envs/env_lofreq.yaml'
    shell:
        """
        lofreq call --call-indels -f {params.ref} -o {output} {input.bam} -q 20 -Q 20 -m 30 -C 3
        """
rule lofreq_call_ignore:
    input:
      bam = rules.mapping_duplicate_samblaster.output.bam#'Sample_{s}/{s}.viral.bam'
    output:
        "Sample_{s}/dedup/lofreq/{s}_lofreq.tsv"
    params:
        ref = config['reference'] 
    conda: 'envs/env_lofreq.yaml'
    shell:
        """
        lofreq call --call-indels -f {params.ref} -o {output} {input.bam} -q 20 -Q 20 -m 30 -C 3
        """
# VarScan
rule varscan:
  input:
      bam = rules.sort_by_position_overlap.output.bamclipoverlap_sorted_bam
  output:
      "Sample_{s}/varscan/{s}_varscan.tsv"
  log:
      "logs/{s}_%s_varscan.log"%datetime.now().strftime('%H_%M_%d_%m_%Y')
  conda: 'envs/env_varscan.yaml'
  #conda: 'envs/env_samtools.yaml'
  params:
      ref = config['reference'],

  shell:
      """
      samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} {input.bam} | varscan pileup2snp --min-reads2 3 --min-coverage 3 --min-avg-qual 20  --variants - > {output}
      """
#previously: samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} {input.bam} | java -jar /path/VarScan.v2.4.6.jar pileup2snp --variants - > {output}
#--p-value 0.01

# VarScan
rule varscan_ignore:
  input:
    bam = rules.mapping_duplicate_samblaster.output.bam
  output:
    "Sample_{s}/dedup/varscan/{s}_varscan.tsv"
  log:
    "logs/{s}_%s_varscan_ignore.log"%datetime.now().strftime('%H_%M_%d_%m_%Y')
  conda: 'envs/env_varscan.yaml'
  #conda: 'envs/env_samtools.yaml'
  params:
    ref = config['reference']
  shell:
    """
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} {input.bam} | varscan pileup2snp --min-reads2 3 --min-coverage 3 --min-avg-qual 20 --variants - > {output}
    """
#--p-value 0.01

# Ivar
rule ivar:
  input:
    bam = rules.sort_by_position_overlap.output.bamclipoverlap_sorted_bam,
    idx = rules.index_overlap.output.idx
  output:
    vcf = 'Sample_{s}/ivar/{s}_ivar.vcf.gz'
  log:
    'Sample_{s}/{s}_vc.log'
  conda: 'envs/env_ivar.yaml'
  params:
    q = 20, #Minimum quality score threshold to count base (Default: 20)
    t = 0,#t=0.03, #Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
    m = 3,
    ref = config['reference'],
    ivar_variants = srcdir("source/ivar_variants_to_vcf.py")
  shell:
    """
    mkdir -p ivar
    cd Sample_{wildcards.s}/ivar
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} ../../{input.bam} | ivar variants -p {wildcards.s}_ivar -q {params.q} -t {params.t} -r {params.ref}
    python {params.ivar_variants} {wildcards.s}_ivar.tsv {wildcards.s}_ivar.vcf > {wildcards.s}.variant.counts.log
    bgzip -c {wildcards.s}_ivar.vcf > {wildcards.s}_ivar.vcf.gz
    """
#    samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} ../virus/{wildcards.s}_position_sorted.bam | ivar variants -p {wildcards.s}_ivar -q {params.q} -t {params.t} -m {params.m} -r {params.ref}

rule ivar_ignore:
  input:
    bam = rules.mapping_duplicate_samblaster.output.bam,
    idx = rules.mapping_duplicate_samblaster.output.idx
  output:
    vcf = 'Sample_{s}/dedup/ivar/{s}_ivar.vcf.gz'
  log:
    'Sample_{s}/{s}_vs_dedup.log'
  conda: 'envs/env_ivar.yaml'
  params:
    q = 20, #Minimum quality score threshold to count base (Default: 20)
    t = 0,#t=0.03, Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
    ref = config['reference'],
    m = 3,
    ivar_variants = srcdir("source/ivar_variants_to_vcf.py")
  shell:
    """
    mkdir -p dedup/ivar
    cd Sample_{wildcards.s}/dedup/ivar
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} ../../../{input.bam} | ivar variants -p {wildcards.s}_ivar -q {params.q} -t {params.t} -m {params.m} -r {params.ref}
    python {params.ivar_variants} {wildcards.s}_ivar.tsv {wildcards.s}_ivar.vcf > {wildcards.s}.variant.counts.log
    bgzip -c {wildcards.s}_ivar.vcf > {wildcards.s}_ivar.vcf.gz
    """


# QC
rule ivar_variant_qc:
  input:
    rules.ivar.output.vcf#'Sample_{s}/virus/{s}.vcf.gz'
  output:
    vcf = 'Sample_{s}/virus/{s}_ivar.qcML'
  conda: 'envs/env_ngs_bits.yaml'
  shell:
    """
	VariantQC -in {input} -out virus/{wildcards.s}.qcML
    """

#_______ CONSENSUS ____________________________________________________#

# Ivar
rule consensus_ivar:
  input:
    #ivar = 'Sample_{s}/.viral.bam'
    ivar  = rules.mapping.output.bam

  output:
    fasta = 'Sample_{s}/consensus/{s}_consensus_ivar.fa'
  conda: 'envs/env_ivar.yaml'
  params:
    min_depth = 10,
    freq_threshold = 0.9
  shell:
    """
    samtools mpileup -aa -A -d 0 -Q 0 {input} | \
        ivar consensus -t {params.freq_threshold} -m {params.min_depth} \
        -p Sample_{wildcards.s}/consensus/{wildcards.s}_consensus_ivar
    """

# Umivar2
rule consensus_umivar:
  input:
    vcf = 'Sample_{s}/cfdna/{s}.cfdna_var.vcf.gz'
  output:
    fasta = 'Sample_{s}/consensus/{s}_umivar.fasta'
  conda: 'envs/env.yaml'
  params:
    ref = config['reference'] 
  shell:
    """
    bcftools norm -f {params.ref} {input.vcf} -Oz -o Sample_{wildcards.s}/cfdna/calls.norm.vcf.gz
    bcftools filter --IndelGap 5 Sample_{wildcards.s}/cfdna/calls.norm.vcf.gz -Oz -o Sample_{wildcards.s}/cfdna/calls.filt.vcf.gz
    tabix Sample_{wildcards.s}/cfdna/calls.filt.vcf.gz
    cat {params.ref} | bcftools consensus Sample_{wildcards.s}/cfdna/calls.filt.vcf.gz > {output.fasta}
    """

rule ivar_variant_table:
  input: 
    consensus=rules.ivar.output,
    vcf='Sample_{s}/ivar/{s}_ivar.vcf.gz'
  output: 'Sample_{s}/ivar/{s}.variants.txt'
  conda:'envs/env.yaml'
  shell:
    """
    echo -e '#chr\tpos\tref\talt\tao\tdp\tfreq' > {output}
    tabix -p vcf -f {input.vcf}
    bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%ALT_FREQ\t%DP]\n' {input.consensus} \
        | awk -v OFS='\t' '$7 = $5/$6' \
        >> {output}
    """

rule ivar_variant_table_ignore:
  input: 
    consensus=rules.ivar_ignore.output,
    vcf='Sample_{s}/dedup/ivar/{s}_ivar.vcf.gz'
  output: 'Sample_{s}/dedup/ivar/{s}.variants.txt'
  conda:'envs/env.yaml'
  shell:
    """
    echo -e '#chr\tpos\tref\talt\tao\tdp\tfreq' > {output}
    tabix -p vcf -f {input.vcf}
    bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%ALT_FREQ\t%DP]\n' {input.consensus} \
        | awk -v OFS='\t' '$7 = $5/$6' \
        >> {output}
    """

#_______ ASSEMBLY ____________________________________________________#

rule assembly:
  input:
    r1=rules.kraken.output.out1,#'Sample_{s}/{s}_viral_1.fastq.gz',
    r2=rules.kraken.output.out2#'Sample_{s}/{s}_viral_2.fastq.gz'
  output:
    dir = 'Sample_{s}/spades/contigs.fasta'
  threads: threads_max
  conda: 'envs/env_spades.yaml'
  shell:"spades.py -1 {input.r1} -2 {input.r2} -t {threads} -o Sample_{wildcards.s}/spades/"

rule quast:
    input:
        contigs = rules.assembly.output, 
        ref = config['reference'] 
    output:directory('Sample_{s}/quast_results')
    conda: 'envs/env_quast.yaml'
    shell:"quast -o {output} -r {input.ref} {input.contigs}"
