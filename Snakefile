import pathlib
import re
from datetime import datetime
import os
import glob

all_sampledirs = [ x.name for x in pathlib.Path().iterdir() \
    if x.is_dir() and x.name.startswith("Sample_") ]
all_sampledirs.sort()
all_sampleids = [ re.sub(r'^Sample_', '', x) for x in all_sampledirs ]
all_sampleids = ["20076a109_02",
"20076a110_02",
"20076a111_01",
"20076a112_02",
"20076a113_02",
"20076a114_01",
"20076a115_01",
"20076a116_01",
"20076a117_02",
"20076a118_02",
"20076a119_01",
"20076a120_01",
"20076a121_01",
"20076a122_01",
"20076a123_01",
"20076a124_01",
"20076a125_01",
"20076a126_01",
"20076a127_01",
"20076a128_01",
"20076a129_01",
"20076a130_01",
"20076a131_01",
"20076a132_01",
"20076a133_01",
"20076a134_01",
"20076a135_01",
"20076a136_01",
"20076a137_01"]
threads_max = 8

ref_path = '%s/ref/MN908947.3.fasta'%os.getcwd()
target_twist = '%s/twist_sars_cov2.bed'%os.getcwd()

def get_fastqs_by_R1(wildcards):
    return glob.glob("%s/Sample_%s/%s*R1*.fastq.gz"%(os.getcwd(),wildcards.s,wildcards.s))
    
def get_fastqs_by_R2(wildcards):
    return glob.glob("%s/Sample_%s/%s*R2*.fastq.gz"%(os.getcwd(),wildcards.s,wildcards.s))

rule all:
  input:
    'ref/MN908947.3.fasta',
    expand('Sample_{s}/virus/{s}_position_sorted.bam', s = all_sampleids), # mapping
    expand('Sample_{s}/virus/{s}_position_sorted.bam.bai', s = all_sampleids), # sort and index
    expand('Sample_{s}/{s}_stats_map.qcML', s = all_sampleids), # mapping_qc
    expand('Sample_{s}/{s}.mosdepth.summary.txt', s = all_sampleids), # coverage depth qc
    expand('Sample_{s}/virus/{s}.vcf.gz', s = all_sampleids), # variant calling
    expand('Sample_{s}/{s}_consensus.fa', s = all_sampleids),  # consensus
    #expand('Sample_{s}/quast_results', s = all_sampleids), # assembly 
    expand('Sample_{s}/{s}.variants.txt', s = all_sampleids), # variants.txt

rule download_ref:
    output:'ref/MN908947.3.fasta'
    conda: 'envs/env_bwa.yaml'
    shell:"""
    mkdir -p ref
    cd ref
    wget https://www.ncbi.nlm.nih.gov/search/api/sequence/MN908947.3/?report=fasta -O MN908947.3.fasta
    bwa index MN908947.3.fasta
    """

#### KRAKEN ####
rule kraken:
    input:
        in1='Sample_{s}/{s}_R1_trimmed.fastq.gz',
        in2='Sample_{s}/{s}_R2_trimmed.fastq.gz'
    output:
        out1='Sample_{s}/{s}_kraken_1.fastq.gz',
        out2='Sample_{s}/{s}_kraken_2.fastq.gz',
    threads: 
        8
    conda: 'envs/env_kraken.yaml'
    params:
        db='/mnt/users/ahkocae1/tools/kraken2/kraken2_human'
    shell:
        "kraken2 --db {params.db} -threads {threads} --unclassified-out Sample_{wildcards.s}/{wildcards.s}_kraken#.fastq.gz --report Sample_{wildcards.s}/{wildcards.s}.kraken2.report.txt --report-zero-counts --paired --gzip-compressed {input.in1} {input.in2}"

#### TRIMMING ####

rule trimming:
    input:
        in1=get_fastqs_by_R1,
        in2=get_fastqs_by_R2
    output:
        out1='Sample_{s}/{s}_R1_trimmed.fastq.gz',
        out2='Sample_{s}/{s}_R2_trimmed.fastq.gz',
        qc= 'Sample_{s}/{s}_stats_fastq.qcML'
    threads: 
        8
    conda: 'envs/env_ngs_bits.yaml'
    params:
        a1='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA',
	a2='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    shell:
        "SeqPurge -in1 {input.in1} -in2 {input.in2} -out1 {output.out1} -out2 {output.out2} -a1 {params.a1} -a2 {params.a2} -qc {output.qc} -threads {threads}"

#### MAPPING ####

rule mapping:
  input: 
     r1='Sample_{s}/{s}_kraken_1.fastq.gz',
     r2='Sample_{s}/{s}_kraken_2.fastq.gz',
     ref= ref_path
  output:
    bam = 'Sample_{s}/virus/{s}_unsorted.bam'
  threads: 8
  conda: 'envs/env_bwa.yaml'
  shell:
    """
    mkdir -p virus
    bwa mem -t {threads} {input.ref} {input.r1} {input.r2} -M | samtools view -bS - > Sample_{wildcards.s}/virus/{wildcards.s}_unsorted.bam
    """
		
rule mapping_qc:
  input:
    bam = 'Sample_{s}/virus/{s}_position_sorted.bam',
    idx = 'Sample_{s}/virus/{s}_position_sorted.bam.bai'
  output:
    qcml = 'Sample_{s}/{s}_stats_map.qcML'
  params:
    target = target_twist
  conda: 'envs/env_ngs_bits.yaml'
  shell:
    """
    MappingQC -out {output.qcml} -roi {params.target} -no_cont -in {input.bam}
    """

rule cov_depth_qc:
  input:
    bam = 'Sample_{s}/virus/{s}_position_sorted.bam',
    idx = 'Sample_{s}/virus/{s}_position_sorted.bam.bai'
  output:
    qcml = 'Sample_{s}/{s}.mosdepth.summary.txt'
  params:
    target = target_twist
  conda: 'envs/env_mosdepth.yaml'
  threads: 2
  shell:
    """
    cd Sample_{wildcards.s}
    mosdepth --by {params.target} -t {threads} {wildcards.s} virus/{wildcards.s}_position_sorted.bam
    """

#### SORT AND INDEX ####

rule sort_by_position:
    input: 
        'Sample_{s}/virus/{s}_unsorted.bam'
    output:
        position_sorted_bam = 'Sample_{s}/virus/{s}_position_sorted.bam'
    conda: 'envs/env_samtools.yaml'
    threads: 8
    shell:"samtools sort -@ {threads} -m 1G -o {output.position_sorted_bam} {input}"	

rule index:
    input: 
        position_sorted_bam = 'Sample_{s}/virus/{s}_position_sorted.bam'
    output:
        position_sorted_idx = 'Sample_{s}/virus/{s}_position_sorted.bam.bai'
    conda: 'envs/env_samtools.yaml'
    threads: 8
    shell:
        "samtools index {input.position_sorted_bam}"	

#### ASSEMBLY ####

rule assembly:
  input:
    r1='Sample_{s}/{s}_kraken_1.fastq.gz',
    r2='Sample_{s}/{s}_kraken_2.fastq.gz'
  output:
    dir = 'Sample_{s}/spades/contigs.fasta'
  threads: 8
  conda: 'envs/env_spades.yaml'
  shell:"spades.py --careful -1 {input.r1} -2 {input.r2} -t {threads} -o {output.dir}"

rule quast:
    input:
        contigs = rules.assembly.output, 
        ref = ref_path
    output:directory('Sample_{s}/quast_results')
    conda: 'envs/env_quast.yaml'
    shell:"quast -o {output} -r {input.ref} {input.contigs}"

#### VARIANT CALLING ####

rule variant_calling:
  input:
    bam = 'Sample_{s}/virus/{s}_position_sorted.bam',
    idx = 'Sample_{s}/virus/{s}_position_sorted.bam.bai',
    ref = ref_path
  output:
    vcf = 'Sample_{s}/virus/{s}.vcf.gz'
  log:
    '{s}_vc.log'
  conda: 'envs/env_ivar.yaml'
  shell:
    """
    cd Sample_{wildcards.s}/virus
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference {input.ref} {wildcards.s}_position_sorted.bam | ivar variants -p {wildcards.s} -q 20 -t 0.03 -r {input.ref}
    python ../../source/ivar_variants_to_vcf.py {wildcards.s}.tsv {wildcards.s}.vcf > {wildcards.s}.variant.counts.log
    bgzip -c {wildcards.s}.vcf > {wildcards.s}.vcf.gz
    """

rule variant_qc:
  input:
    'Sample_{s}/virus/{s}.vcf.gz'
  output:
    vcf = 'Sample_{s}/{s}_variant_calling.qcML'
  conda: 'envs/env_ngs_bits.yaml'
  shell:
    """
	VariantQC -in {input} -out {wildcards.s}.qcML
    """

#### CONSENSUS ####

rule consensus:
  input:
    vcf = 'Sample_{s}/virus/{s}_position_sorted.bam',
  output:
    fasta = 'Sample_{s}/{s}_consensus.fa'
  conda: 'envs/env_ivar.yaml'
  params:
    min_depth = 10,
    freq_threshold = 0.9
  shell:
    """
    samtools mpileup -aa -A -d 0 -Q 0 {input} | \
        ivar consensus -t {params.freq_threshold} -m {params.min_depth} \
        -p Sample_{wildcards.s}/{wildcards.s}_consensus
    """

rule variant_table:
  input: 
    consensus=rules.variant_calling.output,
    vcf='Sample_{s}/virus/{s}.vcf.gz'
  output: 'Sample_{s}/{s}.variants.txt'
  conda:'envs/env.yaml'
  shell:
    """
    echo -e '#chr\tpos\tref\talt\tao\tdp\tfreq' > {output}
    tabix -p vcf -f {input.vcf}
    bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%ALT_FREQ\t%DP]\n' {input.consensus} \
        | awk -v OFS='\t' '$7 = $5/$6' \
        >> {output}
    """
