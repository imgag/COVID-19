# COVID-19

## How to run snakefiles 

- Copy the files and directories from the repository to the project directory where the examples are located.

```
    root
	    |
	    |_ [PROJECT_NAME]
			    |
			    |_ Sample_RXxxxx1
			    |_ Sample_RXxxxx2
			    |_ ....
			    |_ ....
			    |_ envs
			    |_ source
			    |_ ref
			    |_ config.yaml
			    |_ twist_sars_cov2.bed
			    |_ Snakefile_UMI
			    |_ Snakefile_Twist
	    |_ kraken
	            |_ kraken2_human
```  

- Download kraken db  
```
wget https://zenodo.org/record/3738199/files/kraken2_human.tar.gz
```
- Extract it 
```
tar -xvzf kraken2_human.tar.gz
```
- Configure the config.yaml file
```
samples: 'RXxxxx1,RXxxxx2'     # list of samples without "Sample_"
target_twist: '/path/to/twist_sars_cov2.bed' # e.g: /root/[PROJECT_NAME]/twist_sars_cov2.bed
reference: '/path/to/ref/MN908947.3.fasta' # e.g: /root/[PROJECT_NAME]/ref/MN908947.3.fasta
adapter1: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'  # Adapter read 1
adapter2: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'  # Adapter read 2 
krakendb: '/root/kraken/kraken2_human'  # extracted kraken path
threads: 8 # number of threads
```
- running snakefile 
```
cd [PROJECT_NAME]
snakemake -s Snakefile_UMI --configfile config.yaml --cores
or
snakemake -s Snakefile_Twist --configfile config.yaml --cores
```