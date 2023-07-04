# Variant Calling
## Methods from ["Evaluating assembly and variant calling software for strain-resolved analysis of large DNA viruses"](https://academic.oup.com/bib/article/22/3/bbaa123/5868070)
### Variant calling Parameters

> Quality-controlled reads were mapped against the reference genome of the HCMV strains Merlin and AD169 using BWA-MEM with a seed length of 31. HCMV Merlin and AD169 genomes were used as reference genomes, as they were the major strains in all mixtures. The resulting BAM files were deduplicated with the Picard package (http://broadinstitute.github.io/picard/) to remove possible amplification duplicates that may bias the allele frequency of identified variants. To compare the performance of different variant callers, we used LoFreq (parameter: -q 20 -Q 20 -m 20), VarScan2 (—min-avg-qual 20 —P-value 0.01), FreeBayes (—p 1 -m 20 -q 20 -F 0.01 —min-coverage 10), CLC (overall read depth ≥10, average basecall quality ≥20, forward/reverse read balance 0.1–0.9 and variant frequency ≥0.1%), BCFtools (—p 0.01 —ploidy 1 -mv -Ob) and GATK HaplotypeCaller (—min-base-quality-score 20 -ploidy 1) to identify variants. The variants from the difference between genomes detected by MUMmer were considered as positive variants. Based on this standard, precision, recall and F1-score were computed to evaluate those callers. The pairwise genome differences of 30 E. coli or 30 HIV genomes were determined by MUMmer as well. To evaluate the performance of different callers for SNP and InDel prediction, the command vcfeval in RTG-tools [60] was used to generate recall-precision curves based on the Phred scaled ‘QUAL’ score field (—squash-ploidy -f QUAL —sample ALT)."

- https://backoffice.biblio.ugent.be/download/01GY7GK7NPE15KKNYZNP8K8XYC/01GY7GQHFYX9WY4GEE0Z874QK0
page 169
> " LoFreq v2.1.3.1 package, setting
the strand bias threshold for reporting a variant to the maximum allowed value by using the
option “--sb-thresh 2147483647” to allow highly strand-biased variants to be retained, to
account for the non-random distribution of reads due to the design of the amplification panel."

1. umivar parameters
```
tumor vs normal file?
-tbam  tumor bam file
-ac 4:
  -ac AC, --ac AC       Minimum number of reads supporting a variant/not same as depth of coverage, is there any depth of coverage????
-ns -1
  -ns NUM_SITES, --num_sites NUM_SITES
                        Number of sites to be analysed
-sb 0
  -sb {0,1}, --strand_bias {0,1}
                        Fisher strand bias filter. Default [0]
  -mq MQ, --mq MQ       Minimum mapping quality
  -bq BQ, --bq BQ       Minimum base quality

usage: umiVar - variant calling with unique molecular barcodes
       [-h] -tbam TBAM [-nbam NBAM] -r REF [-b BED] [-m MONITORING]
       [-o OUT_FOLDER] [-p PARAM] [-mq MQ] [-bq BQ] [-d DIST] [-ac AC]
       [-af AF] [-ns NUM_SITES] [-sb {0,1}] [-t TEMP_DIR] [-kt]
       [-crb CUSTOM_RSCRIPT_BINARY]

optional arguments:
  -h, --help            show this help message and exit
  -tbam TBAM, --tbam TBAM
                        Tumor bam file
  -nbam NBAM, --nbam NBAM
                        Normal bam file
  -r REF, --ref REF     Reference genome - fasta
  -b BED, --bed BED     Bed file of the targeted regions. O-based
  -m MONITORING, --monitoring MONITORING
                        VCF file with genomic positions for monitoring or
                        sample-IDing.
  -o OUT_FOLDER, --out_folder OUT_FOLDER
                        Output folder. Will be created if not existing
  -p PARAM, --param PARAM
                        Beta-binomial parameters table
  -mq MQ, --mq MQ       Minimum mapping quality
  -bq BQ, --bq BQ       Minimum base quality
  -d DIST, --dist DIST  Minimum distance between variants
  -ac AC, --ac AC       Minimum number of reads supporting a variant
  -af AF, --af AF       Minimum fraction of reads supporting a variant
  -ns NUM_SITES, --num_sites NUM_SITES
                        Number of sites to be analysed
  -sb {0,1}, --strand_bias {0,1}
                        Fisher strand bias filter. Default [0]
  -t TEMP_DIR, --temp_dir TEMP_DIR
                        Temporary directory
  -kt, --keep_temp      Don't delete temporary directory

```
2. varscan parameters
```
> VarScan2 (—min-avg-qual 20 —P-value 0.01)

samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} {input.bam} | varscan pileup2snp --variants - > {output}
samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} {input.bam} | varscan pileup2snp --min-reads 4 --min-coverage 4 --min-avg-qual 20 --p-value 0.01 --variants - > {output}

varscan pileup2snp -h         
Warning: No p-value threshold provided, so p-values will not be calculated                                                               
Min coverage:   8                                                                                                                        
Min reads2:     2                                                                                                                        
Min var freq:   0.01                                                                                                                     
Min avg qual:   15                                                                                                                       
P-value thresh: 0.99                                                                                                                     
USAGE: java -jar VarScan.jar pileup2cns [pileup file] OPTIONS                                                                            
        pileup file - The SAMtools pileup file                                                                                           
                                                                                                                                
        OPTIONS:                                                                                                                         
        --min-coverage  Minimum read depth at a position to make a call [8]                                                              
        --min-reads2    Minimum supporting reads at a position to call variants [2]                                                      
        --min-avg-qual  Minimum base quality at a position to count a read [15]                                                          
        --min-var-freq  Minimum variant allele frequency threshold [0.01]                                                                
        --min-freq-for-hom      Minimum frequency to call homozygote [0.75]                                                              
        --p-value       Default p-value threshold for calling variants [99e-02]                                                          
        --variants      Report only variant (SNP/indel) positions [0] 
```
3. lofreq
call variants from BAM file . LoFreq (parameter: -q 20 -Q 20 -m 20)
pipeline parameter:
lofreq call --call-indels -f /mnt/storage2/users/ahcepev1/pipelines/COVID-19/ref/MN908947.3.fasta -o Sample_21014a009_01/lofreq/21014a009_01_lofreq.tsv Sample_21014a009_01/dedup/21014a009_01_bamclipoverlap_sorted.bam

```
Options:
- Reference:
       -f | --ref FILE              Indexed reference fasta file (gzip supported) [null]
- Output:
       -o | --out FILE              Vcf output file [- = stdout]
- Base-call quality:
       -q | --min-bq INT            Skip any base with baseQ smaller than INT [6]
       -Q | --min-alt-bq INT        Skip alternate bases with baseQ smaller than INT [6]
- Mapping quality:
       -m | --min-mq INT            Skip reads with mapping quality smaller than INT [0]
- P-values:
       -a | --sig                   P-Value cutoff / significance level [0.010000]
       -b | --bonf                  Bonferroni factor. 'dynamic' (increase per actually performed test) or INT ['dynamic']
- Indels:1
            --call-indels           Enable indel calls (note: preprocess your file to include indel alignment qualities!)
            --only-indels           Only call indels; no SNVs
- Misc.:
       -C | --min-cov INT           Test only positions having at least this coverage [1]
                                    (note: without --no-default-filter default filters (incl. coverage) kick in after predictions are done)
       -d | --max-depth INT         Cap coverage at this depth [1000000]

```

4. ivar
There are two parameters that can be set for variant calling using iVar - minimum quality(Default: 20) and minimum frequency(Default: 0.03). Minimum quality is the minimum quality for a base to be used in frequency calculations at a given position. Minimum frequency is the minimum frequency required for a SNV or indel to be reported.

"For Samtools mpileup mostly default parameters are used. The only non-default parameters used are do not discard anomalous read pairs (-A), disable per-base alignment quality (-B), skip bases with base quality smaller than 0 (-Q 0), and especially the maximum per-file coverage -d 1000000. The coverage of the pool amplicons can vary more than 100x in tiled amplicon approaches. Therefore, the value of the -d parameter should be at least 5-10x greater than the Avg. Coverage (Unassembled). For consensus calling by default the minimum quality score threshold (-q 20), minimum coverage to call consensus (-m 10), and minimum frequency threshold (-t 0.7) parameter values are applied. For stricter consensus calling, e.g., a called base must make up at least 90% presence at a position, the latter parameter must be changed to -t 0.9. For further details and other parameters it is referred to the iVar manual. All by SeqSphere+" --https://www.ridom.de/u/SARS-CoV-2_Analysis_Quick_Start.html
```
params:
    q = 20, #Minimum quality score threshold to count base (Default: 20)
    t = 0.03, #Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
    ref = config['reference'],
    ivar_variants = srcdir("source/ivar_variants_to_vcf.py")
  shell:
    """
    mkdir -p ivar
    cd Sample_{wildcards.s}/ivar
    samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} ../../{input.bam} | ivar variants -p {wildcards.s}_ivar -q {params.q} -t {params.t} -r {params.ref}
    python {params.ivar_variants} {wildcards.s}_ivar.tsv {wildcards.s}_ivar.vcf > {wildcards.s}.variant.counts.log
    bgzip -c {wildcards.s}_ivar.vcf > {wildcards.s}_ivar.vcf.gz
```
-m   4 Minimum read depth to call variants (Default: 0)

ivar variants
Usage: samtools mpileup -aa -A -d 0 -B -Q 0 --reference [<reference-fasta] <input.bam> | ivar variants -p <prefix> [-q <min-quality>] [-t <min-frequency-threshold>] [-m <minimum depth>] [-r <reference-fasta>] [-g GFF file]
Note : samtools mpileup output must be piped into ivar variants

Input Options    Description
           -q    Minimum quality score threshold to count base (Default: 20)
           -t    Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)
           -m    Minimum read depth to call variants (Default: 0)
           -r    Reference file used for alignment. This is used to translate the nucleotide sequences and identify intra host single nucleotide variants
           -g    A GFF file in the GFF3 format can be supplied to specify coordinates of open reading frames (ORFs). In absence of GFF file, amino acid translation will not be done.

Output Options   Description
           -p    (Required) Prefix for the output tsv variant file
```
