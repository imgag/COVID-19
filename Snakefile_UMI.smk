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
    fastqs = glob.glob("Sample_%s/%s*R1*.fastq.gz"%(wildcards.s,wildcards.s))
    return sorted(fastqs,key=lambda x: x)
    
def get_fastqs_by_R2(wildcards):
    fastqs = glob.glob("Sample_%s/%s*R2*.fastq.gz"%(wildcards.s,wildcards.s))
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
    expand('Sample_{s}/virus/{s}_stats_map.qcML', s = all_sampleids), # mapping_qc
    expand('Sample_{s}/coverage/{s}.mosdepth.summary.txt', s = all_sampleids), # coverage depth qc
    expand('Sample_{s}/ivar/{s}_ivar.vcf.gz', s = all_sampleids), # variant calling
    expand('Sample_{s}/consensus/{s}_consensus_ivar.fa', s = all_sampleids),  # consensus ivar
    #expand('Sample_{s}/consensus/{s}_umivar.fasta', s = all_sampleids),  # consensus umivar
    #expand('Sample_{s}/quast_results', s = all_sampleids), # assembly 
    expand('Sample_{s}/ivar/{s}.variants.txt', s = all_sampleids), # variants.txt
    expand('Sample_{s}/umivar2/{s}.viral_hq.vcf', s = all_sampleids), #umiVar2
    expand('Sample_{s}/lofreq/{s}_lofreq.tsv', s = all_sampleids), # Lofreq
    expand('Sample_{s}/varscan/{s}_varscan.tsv', s = all_sampleids) # varscan



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
        out1=temp('Sample_{s}/{s}_R1_barcode_added.fastq.gz'),
        out2=temp('Sample_{s}/{s}_R2_barcode_added.fastq.gz')
    conda: 'envs/env_ngs_bits.yaml'
    shell:
        """
        FastqAddBarcode -in1 {input.in1} -in2 {input.in2} -in_barcode {input.barcode} -out1 {output.out1} -out2 {output.out2}
        """

# Read trimming and ReadQC
rule trimming:
    input:
        in1='Sample_{s}/{s}_R1_barcode_added.fastq.gz',
        in2='Sample_{s}/{s}_R2_barcode_added.fastq.gz'
    output:
        out1=temp('Sample_{s}/{s}_R1_trimmed.fastq.gz'),
        out2=temp('Sample_{s}/{s}_R2_trimmed.fastq.gz'),
        qc= 'Sample_{s}/{s}_stats_fastq.qcML'
    threads: 
        threads_max
    conda: 'envs/env_ngs_bits.yaml'
    params:
        a1=config['adapter1'],
	a2=config['adapter2']
    shell:
        """
        SeqPurge -progress 1000 -in1 {input.in1} -in2 {input.in2} -out1 {output.out1} -out2 {output.out2} -a1 {params.a1} -a2 {params.a2} -qc {output.qc} -threads {threads}
        """

# Removal of human host reads with Kraken 
rule kraken:
    input:
        in1='Sample_{s}/temp/{s}_R1_trimmed.fastq.gz',
        in2='Sample_{s}/temp/{s}_R2_trimmed.fastq.gz'
    output:
        out1='Sample_{s}/{s}_viral_R1.fastq.gz',
        out2='Sample_{s}/{s}_viral_R2.fastq.gz',
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
        --unclassified-out Sample_{wildcards.s}/{wildcards.s}_viral_R#.fastq  \
        --report Sample_{wildcards.s}/{wildcards.s}.kraken2.report.txt  \
        --report-zero-counts  \
        --paired --gzip-compressed \
        {input.in1} {input.in2} > /dev/null
       gzip Sample_{wildcards.s}/{wildcards.s}_viral_R1.fastq
       gzip Sample_{wildcards.s}/{wildcards.s}_viral_R2.fastq
       """


#_______ MAPPING ___________________________________________________________#

rule mapping:
  input: 
     r1='Sample_{s}/{s}_viral_R1.fastq.gz',
     r2='Sample_{s}/{s}_viral_R2.fastq.gz',
     index_files = rules.index_reference.output
  output:
    bam = 'Sample_{s}/{s}.viral.bam'
  threads: threads_max
  params: 
    ref = config['reference']
  conda: 'envs/env_bwa.yaml'
  shell:
    """
    bwa mem -t {threads} {params.ref} {input.r1} {input.r2} -M | \
      samblaster -M | \
      samtools view -bS - | \
      samtools sort -@4 -m 1G -o {output} -
    samtools index {output}
    """
		
rule mapping_qc:
  input:
    bam = 'Sample_{s}/{s}.viral.bam'
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
    bam = 'Sample_{s}/{s}.viral.bam'
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

#_______ DEDUPLICATE ________________________________________________________#

rule barcode_correction:
  input:
    bam = 'Sample_{s}/{s}.viral.bam'
  output: 
    corrected_bam = 'Sample_{s}/dedup/{s}.viral.corrected.bam'
  conda:
     "envs/env_umivar.yaml"
  log:
    "logs/{s}_barcode_correction.log"    
  params:
    barcode_correction = '%s/barcode_correction.py'%config['umiVar']
  shell:
    """
    {params.barcode_correction} --infile {input.bam} --outfile {output.corrected_bam} > {log}
    """

rule sort_by_position_after_correction:
    input: 
        'Sample_{s}/dedup/{s}_corrected.bam'
    output:
        corrected_sorted_bam = 'Sample_{s}/dedup/{s}_corrected_sorted.bam'
    conda: 'envs/env_samtools.yaml'
    threads: threads_max
    shell:"samtools sort -@ {threads} -m 1G -o {output.corrected_sorted_bam} {input}"	

rule index_after_correction:
    input: 
        corrected_sorted_bam = 'Sample_{s}/dedup/{s}_corrected_sorted.bam'
    output:
        'Sample_{s}/dedup/{s}_corrected_sorted.bam.bai'
    conda: 'envs/env_samtools.yaml'
    shell:
        "samtools index {input.corrected_sorted_bam}"	

rule bamclipoverlap:
    input:
       corrected_sorted_bam = 'Sample_{s}/dedup/{s}_corrected_sorted.bam',
       corrected_sorted_idx = 'Sample_{s}/dedup/{s}_corrected_sorted.bam.bai'
    output:
       overlapped_bam = 'Sample_{s}/dedup/{s}_bamclipoverlap.bam'
    conda: 'envs/env_ngs_bits.yaml'
    shell: 
       """
	   BamClipOverlap -in {input.corrected_sorted_bam} -out {output.overlapped_bam} -overlap_mismatch_basen
       """

rule sort_by_position_overlap:
    input: 
        'Sample_{s}/dedup/{s}_bamclipoverlap.bam'
    output:
        bamclipoverlap_sorted_bam = 'Sample_{s}/dedup/{s}_bamclipoverlap_sorted.bam'
    conda: 'envs/env_samtools.yaml'
    threads: threads_max
    shell:"samtools sort -@ {threads} -m 1G -o {output.bamclipoverlap_sorted_bam} {input}"	

rule index_overlap:
    input: 
        bamclipoverlap_sorted_bam = 'Sample_{s}/dedup/{s}_bamclipoverlap_sorted.bam'
    output:
        bamclipoverlap_sorted_idx = 'Sample_{s}/dedup/{s}_bamclipoverlap_sorted.bam.bai'
    conda: 'envs/env_samtools.yaml'
    shell:
        "samtools index {input.bamclipoverlap_sorted_bam}"	

#_______ VARIANT CALLING ________________________________________________________#

# Umivar2
rule umivar2:
  input: 
    bam = 'Sample_{s}/{s}.viral.bam'
  output: 'Sample_{s}/umivar2/{s}.viral_hq.vcf'
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
    kt = "-kt " if 'kt' in config['umivar'] else ""

  shell:
    """
    {params.umiVar} \
        -tbam {input.bam} \
        -b {params.target} \
        -r {params.ref} \
        -o /mnt/users/ahgrosc1/projects/covid_umi/samples/Sample_{wildcards.s}/umivar2 \
        {params.ac} {params.af} {params.ns} {params.sb} {params.kt} 
    """

# LoFreq
rule lofreq_call:
    input:
        bam = 'Sample_{s}/{s}.viral.bam'
    output:
        "Sample_{s}/lofreq/{s}_lofreq.tsv"
    params:
        ref = config['reference'] 
    conda: 'envs/env_lofreq.yaml'
    shell:
        """
        lofreq call --call-indels -f {params.ref} -o {output} {input.bam}
        """

# VarScan
rule varscan:
  input:
    bam = 'Sample_{s}/{s}.viral.bam'
  output:
    "Sample_{s}/varscan/{s}_varscan.tsv"
  log:
    "logs/{s}_%s_varscan.log"%datetime.now().strftime('%H_%M_%d_%m_%Y')
  conda: 'envs/env_samtools.yaml'
  params:
    ref = config['reference']
  shell:
    """
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} {input.bam} | java -jar /mnt/share/opt/VarScan.v2.4.4/VarScan.v2.4.4.jar pileup2snp --variants - > {output}
    """
		
# Ivar
rule variant_calling:
  input:
    bam = 'Sample_{s}/{s}.viral.bam',
    idx = 'Sample_{s}/{s}.viral.bam.bai',
  output:
    vcf = 'Sample_{s}/ivar/{s}_ivar.vcf.gz'
  log:
    'Sample_{s}/{s}_vc.log'
  conda: 'envs/env_ivar.yaml'
  params:
    q = 20, #Minimum quality score threshold to count base (Default: 20)
    t = 0.03, #Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
    ref = config['reference'],
    ivar_variants = srcdir("source/ivar_variants_to_vcf.py")
  shell:
    """
    mkdir -p ivar
    cd Sample_{wildcards.s}/ivar
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} ../virus/{wildcards.s}_position_sorted.bam | ivar variants -p {wildcards.s}_ivar -q {params.q} -t {params.t} -r {params.ref}
    python {params.ivar_variants} {wildcards.s}_ivar.tsv {wildcards.s}_ivar.vcf > {wildcards.s}.variant.counts.log
    bgzip -c {wildcards.s}_ivar.vcf > {wildcards.s}_ivar.vcf.gz
    """

# QC
rule variant_qc:
  input:
    'Sample_{s}/virus/{s}.vcf.gz'
  output:
    vcf = 'Sample_{s}/virus/{s}_variant_calling.qcML'
  conda: 'envs/env_ngs_bits.yaml'
  shell:
    """
	VariantQC -in {input} -out virus/{wildcards.s}.qcML
    """

#_______ CONSENSUS ____________________________________________________#

# Ivar
rule consensus_ivar:
  input:
    ivar = 'Sample_{s}/.viral..bam'
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

rule variant_table:
  input: 
    consensus=rules.variant_calling.output,
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
	
#_______ ASSEMBLY ____________________________________________________#

rule assembly:
  input:
    r1='Sample_{s}/{s}_viral_R1.fastq.gz',
    r2='Sample_{s}/{s}_viral_R2.fastq.gz'
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