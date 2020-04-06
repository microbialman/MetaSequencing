#Snakefile to run full MetaSequencing pipeline 
import os

#get the default config files
#this can be overwritten at runtime using --configfile
configfile: workflow.basedir+"/config/global.yaml"

#patterns to define input files
spat="[a-zA-Z0-9_]*"
epat="(fasta$|fasta.gz$|fasta.1.gz$|fasta.1$|fna$|fna.gz$|fna.1.gz$|fna.1$|fa$|fa.gz$|fa.1.gz$|fa.1$|fastq$|fastq.gz$|fastq.1.gz$|fastq.1$)"
inglob="./{sample, "+spat+"}.{ext, "+epat+"}"
wildcard_constraints:
    sample=spat,
    ext=epat

#get the input files from directory and generate dic mapping file to extension
samples,exts=glob_wildcards(inglob)
fileexts={samples[x]:exts[x] for x in range(len(samples))}

#get directory
cwd=os.getcwd()

#the target files, these will define which steps are run
rule all:
    input:
        "Filtering/filter_report.html",
        "Assembly/assemble_report.html",
        expand("Annotation/functional_annotations.dir/emapper_chunks/{sample}.chunk.log",sample=samples)
        
#rule files for each stage
include: "workflow/rules/Filter.smk"
include: "workflow/rules/Assemble.smk"
include: "workflow/rules/Annotate.smk"
