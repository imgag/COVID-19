- parameters used

```
[Paper](https://academic.oup.com/bib/article/22/3/bbaa123/5868070)

"Variant calling

Quality-controlled reads were mapped against the reference genome of the HCMV strains Merlin and AD169 using BWA-MEM with a seed length of 31. HCMV Merlin and AD169 genomes were used as reference genomes, as they were the major strains in all mixtures. The resulting BAM files were deduplicated with the Picard package (http://broadinstitute.github.io/picard/) to remove possible amplification duplicates that may bias the allele frequency of identified variants. To compare the performance of different variant callers, we used LoFreq (parameter: -q 20 -Q 20 -m 20), VarScan2 (—min-avg-qual 20 —P-value 0.01), FreeBayes (—p 1 -m 20 -q 20 -F 0.01 —min-coverage 10), CLC (overall read depth ≥10, average basecall quality ≥20, forward/reverse read balance 0.1–0.9 and variant frequency ≥0.1%), BCFtools (—p 0.01 —ploidy 1 -mv -Ob) and GATK HaplotypeCaller (—min-base-quality-score 20 -ploidy 1) to identify variants. The variants from the difference between genomes detected by MUMmer were considered as positive variants. Based on this standard, precision, recall and F1-score were computed to evaluate those callers. The pairwise genome differences of 30 E. coli or 30 HIV genomes were determined by MUMmer as well. To evaluate the performance of different callers for SNP and InDel prediction, the command vcfeval in RTG-tools [60] was used to generate recall-precision curves based on the Phred scaled ‘QUAL’ score field (—squash-ploidy -f QUAL —sample ALT)."


https://backoffice.biblio.ugent.be/download/01GY7GK7NPE15KKNYZNP8K8XYC/01GY7GQHFYX9WY4GEE0Z874QK0
page 169
" LoFreq v2.1.3.1 package, setting
the strand bias threshold for reporting a variant to the maximum allowed value by using the
option “--sb-thresh 2147483647” to allow highly strand-biased variants to be retained, to
account for the non-random distribution of reads due to the design of the amplification panel."
```
- umivar parameters
```
-tbam  tumor bam file
-ac 4:
  -ac AC, --ac AC       Minimum number of reads supporting a variant

-ns -1
  -ns NUM_SITES, --num_sites NUM_SITES
                        Number of sites to be analysed
-sb 0
  -sb {0,1}, --strand_bias {0,1}
                        Fisher strand bias filter. Default [0]

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
  -crb CUSTOM_RSCRIPT_BINARY, --custom_rscript_binary CUSTOM_RSCRIPT_BINARY
                        Path to custom Rscript binary. [Default: 'Rscript']
```
- varscan parameters
```
samtools mpileup -aa -A -d 0 -B -Q 0 --reference {params.ref} {input.bam} | varscan pileup2snp --variants - > {output}
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
- lofreq

```
lofreq call --call-indels -f /mnt/storage2/users/ahcepev1/pipelines/COVID-19/ref/MN908947.3.fasta -o Sample_21014a009_01/lofreq/21014a009_01_lofreq.tsv Sample_21014a009_01/dedup/21014a009_01_bamclipoverlap_sorted.bam

lofreq call: call variants from BAM file

Usage: lofreq call [options] in.bam

Options:
- Reference:
       -f | --ref FILE              Indexed reference fasta file (gzip supported) [null]
- Output:
       -o | --out FILE              Vcf output file [- = stdout]
- Regions:
       -r | --region STR            Limit calls to this region (chrom:start-end) [null]
       -l | --bed FILE              List of positions (chr pos) or regions (BED) [null]
- Base-call quality:
       -q | --min-bq INT            Skip any base with baseQ smaller than INT [6]
       -Q | --min-alt-bq INT        Skip alternate bases with baseQ smaller than INT [6]
       -R | --def-alt-bq INT        Overwrite baseQs of alternate bases (that passed bq filter) with this value (-1: use median ref-bq; 0: keep) [0]
       -j | --min-jq INT            Skip any base with joinedQ smaller than INT [0]
       -J | --min-alt-jq INT        Skip alternate bases with joinedQ smaller than INT [0]
       -K | --def-alt-jq INT        Overwrite joinedQs of alternate bases (that passed jq filter) with this value (-1: use median ref-bq; 0: keep) [0]
- Base-alignment (BAQ) and indel-aligment (IDAQ) qualities:
       -B | --no-baq                Disable use of base-alignment quality (BAQ)
       -A | --no-idaq               Don't use IDAQ values (NOT recommended under ANY circumstances other than debugging)
       -D | --del-baq               Delete pre-existing BAQ values, i.e. compute even if already present in BAM
       -e | --no-ext-baq            Use 'normal' BAQ (samtools default) instead of extended BAQ (both computed on the fly if not already present in lb tag)
- Mapping quality:
       -m | --min-mq INT            Skip reads with mapping quality smaller than INT [0]
       -M | --max-mq INT            Cap mapping quality at INT [255]
       -N | --no-mq                 Don't merge mapping quality in LoFreq's model
- Indels:
            --call-indels           Enable indel calls (note: preprocess your file to include indel alignment qualities!)
            --only-indels           Only call indels; no SNVs
- Source quality:
       -s | --src-qual              Enable computation of source quality
       -S | --ign-vcf FILE          Ignore variants in this vcf file for source quality computation. Multiple files can be given separated by commas
       -T | --def-nm-q INT          If >= 0, then replace non-match base qualities with this default value [-1]
- P-values:
       -a | --sig                   P-Value cutoff / significance level [0.010000]
       -b | --bonf                  Bonferroni factor. 'dynamic' (increase per actually performed test) or INT ['dynamic']
- Misc.:
       -C | --min-cov INT           Test only positions having at least this coverage [1]
                                    (note: without --no-default-filter default filters (incl. coverage) kick in after predictions are done)
       -d | --max-depth INT         Cap coverage at this depth [1000000]
            --illumina-1.3          Assume the quality is Illumina-1.3-1.7/ASCII+64 encoded
            --use-orphan            Count anomalous read pairs (i.e. where mate is not aligned properly)
            --plp-summary-only      No variant calling. Just output pileup summary per column
            --no-default-filter     Don't run default 'lofreq filter' automatically after calling variants
            --force-overwrite       Overwrite any existing output
            --verbose               Be verbose
```

- ivar

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
