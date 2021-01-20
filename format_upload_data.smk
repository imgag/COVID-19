import datetime
import pandas as pd
configfile: "config.yaml"

li = [pd.read_csv(x) for x in config['data_files']]
dt = pd.concat(li, axis = 0, ignore_index=True)
dt['sampling_date']= pd.to_datetime(dt['sampling_date'], format = "%d.%m.%Y")

print(dt)

def get_fasta(wc):
    f = config['project_prefix'] + "/" + dt[dt.external_sample==wc.s].project.item() + "/" + config['fasta_filename_template']
    f = f.replace("{sample}", dt[dt.external_sample==wc.s].internal_sample.item())
    return(f)

# Return Strain name for wildcards and config file
def get_name(wc):
    year = dt[dt.external_sample==wc.s].sampling_date.item().year
    strain = config['nextstrain_metadata']['strain'].replace('{year}', str(year))
    strain = strain.replace('{sample}', wc.s)
    return(strain)

def get_meta_gisaid(wc):
    date = dt[dt.external_sample==wc.s].sampling_date.item().date
    conv = lambda i : i or '' 
    s = ",".join([str(conv(x)) for x in config['gisaid_metadata'].values()])
    s = s.replace('{sample}', wc.s)
    s = s.replace('{batch}', str(dt[dt.external_sample==wc.s].batch))
    s = s.replace('{year}', str(dt[dt.external_sample==wc.s].sampling_date.item().year))
    s = s.replace('{month}', str(dt[dt.external_sample==wc.s].sampling_date.item().month))
    s = s.replace('{day}', str(dt[dt.external_sample==wc.s].sampling_date.item().day))
    s = s.replace('{today}', str(datetime.date.today()))
    #print(s)
    return(s)

def get_header_gisaid(wc):
    conv = lambda i : i or 'q' 
    s = ",".join([str(conv(x)) for x in config['gisaid_metadata'].keys()])
    return(s)

def get_meta_nextstrain(wc):
    date = dt[dt.external_sample==wc.s].sampling_date.item().date
    conv = lambda i : i or '' 
    s = "\t".join([str(conv(x)) for x in config['nextstrain_metadata'].values()])
    s = s.replace('{sample}', wc.s)
    s = s.replace('{year}', str(dt[dt.external_sample==wc.s].sampling_date.item().year))
    s = s.replace('{month}', str(dt[dt.external_sample==wc.s].sampling_date.item().month))
    s = s.replace('{day}', str(dt[dt.external_sample==wc.s].sampling_date.item().day))
    s = s.replace('{today}', str(datetime.date.today()))
    #print(s)
    return(s)

def get_header_nexstrain(wc):
    conv = lambda i : i or '' 
    s = "\t".join([str(conv(x)) for x in config['nextstrain_metadata'].keys()])
    return(s)

rule all:
    input:
        expand(config['gisaid_metadata']['fn'], batch = dt.batch.unique()),
        expand(config['filename_gisaid_metadata'], batch = dt.batch.unique()),
        expand(config['filename_nextstrain_metadata'], batch = dt.batch.unique())

rule format_fasta:
    input:
        get_fasta
    output:
        "{batch}/fasta/{s}.fasta"
    params:
        sample = get_name
    shell:
        """
        echo ">{params.sample}" > {output}
        tail -n +2 {input} >> {output}
        """

rule merge_fasta:
    input:
        expand("{{batch}}/fasta/{s}.fasta", s = dt['external_sample'])
    output:
        config['gisaid_metadata']['fn']
    shell:
        """
        cat {input} > {output}
        """

rule meta_gisaid:
    output:
        "{batch}/meta_gisaid/{s}.csv"
    params:
        meta = get_meta_gisaid
    shell:
        """
        echo "{params.meta}" > {output}
        """
        
rule merge_meta_gisaid:
    input:
        expand("{{batch}}/meta_gisaid/{s}.csv", s = dt['external_sample'])
    output:
        config['filename_gisaid_metadata']
    params:
        header = get_header_gisaid
    shell:
        """
        echo "{params.header}" > {output}
        cat {input} >> {output}
        """
rule meta_nextstrain:
    output:
        "{batch}/meta_nextstrain/{s}.tsv"
    params:
        meta = get_meta_nextstrain
    shell:
        """
       echo "{params.meta}" > {output}
        """

rule merge_meta_nextstrain:
    input:
        expand("{{batch}}/meta_nextstrain/{s}.tsv", s = dt['external_sample'])
    output:
        config['filename_nextstrain_metadata']
    params:
        header = get_header_nexstrain
    shell:
        """
        echo "{params.header}" > {output}
        cat {input} >> {output}
        """