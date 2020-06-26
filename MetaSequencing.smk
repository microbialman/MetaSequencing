#Snakefile to run full MetaSequencing pipeline 
import os

#get the default config files
#this can be overwritten at runtime using --configfile
configfile: workflow.basedir+"/config/global.yaml"
            
#patterns to define input files
spat="[a-zA-Z0-9_]*"
epat="(fasta$|fasta.gz$|1.fasta.gz$|1.fasta$|fna$|fna.gz$|1.fna.gz$|1.fna$|fa$|fa.gz$|1.fa.gz$|1.fa$|fastq$|fastq.gz$|1.fastq.gz$|1.fastq$|fq$|fq.gz$|1.fq.gz$|1.fq$|contigs.fa.gz$)"
inglob="./{sample, "+spat+"}.{ext, "+epat+"}"
wildcard_constraints:
    sample=spat,
    ext=epat

#get the input files from directory and generate dic mapping file to extension
samples,exts=glob_wildcards(inglob)
fileexts={samples[x]:exts[x] for x in range(len(samples))}

#get directory
cwd=os.getcwd()

#function to make the target files base on stages to be run
def infun(wc):    
    if config["Global"]["run"] == "filter":
        return("Filtering/filter_report.html")
    elif config["Global"]["run"] == "assemble":
        return("Assembly/assemble_report.html")
    elif config["Global"]["run"] == "annotate":
        return("Annotation/annotate_report.html")
    elif config["Global"]["run"] == "enumerate":
        return("Enumeration/enumerate_report.html")
    else:
        return("Filtering/filter_report.html",
        "Assembly/assemble_report.html",
        "Annotation/annotate_report.html",
        "Enumeration/enumerate_report.html")
    
#the target files, these will define which steps are run
rule all:
    input:
        infun

        
#rule files for each stage
include: "workflow/rules/Filter.smk"
include: "workflow/rules/Assemble.smk"
include: "workflow/rules/Annotate.smk"
include: "workflow/rules/Enumerate.smk"
