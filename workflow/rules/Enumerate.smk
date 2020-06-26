#Rules to enumerate annotations from reads and generate gmt files for set ernichment analyses
from workflow.functions import Hisat2 as H
from workflow.functions import FeatureCounts as F

configfile: workflow.basedir+"/config/enumerate.yaml"

#functions to allow switching between enumeration  only or following whole run
def indexin(wc):
    if config["Global"]["run"]=="enumerate":
        return(config["Enumerate"]["General"]["contig_dir"]+"/{sample}"+config["Enumerate"]["General"]["contig_ext"])
    else:
        return("Assembly/contigs.dir/{sample}.contigs.fa.gz")
    
    
#build hisat2 indexes for the contig files
rule makehisatindex:
    input:
        indexin
    output:
        "Enumeration/contig_indices.dir/{sample}/{sample}.1.ht2l"
    threads:
        int(config["Enumerate"]["HISAT2"]["build_threads"])
    resources:
        mem_mb=int(config["Enumerate"]["HISAT2"]["build_memory"])
    run:
        command=H.Hisat2Build(input[0],output[0],config["Enumerate"]).build()
        shell(command)

#function toggling whole pipeline input and enumeration only
def readsin(wc):
    if config["Global"]["run"]=="enumerate":
        return("{{sample}}.{}".format(fileexts[wc.sample]))
    else:
        return("Filtering/genome_filter_out.dir/unmapped/{{sample}}.{}".format(fileexts[wc.sample]))
    
        
#map the reads against the hisat indices
config["Enumerate"]["HISAT2"]["delete_mapped"]="false"
config["Enumerate"]["HISAT2"]["convert_to_bam"]="true"
rule mapreads:
    input:
        readsin,
        "Enumeration/contig_indices.dir/{sample}/{sample}.1.ht2l"
    output:
        temp("Enumeration/mapped_reads.dir/mapped/{sample}.bam")
    threads:
        int(config["Enumerate"]["HISAT2"]["threads"])
    resources:
        mem_mb=int(config["Enumerate"]["HISAT2"]["memory"])
    run:
        index=input[1].replace(".1.ht2l","")
        command=H.Hisat2(input[0],output[0],config["Enumerate"],index).build()
        shell(command)

#rule to toggle input from whole run or just running enumerate
def gtfin(wc):
    if config["Global"]["run"]=="enumerate":
        return(config["Enumerate"]["General"]["gtf_dir"]+"/{sample}.orf_annotations.short.gtf.gz")
    else:
        return("Annotation/combined_annotations.dir/{sample}.orf_annotations.short.gtf.gz")

#count the orfs from the mapped reads
rule runfeaturecounts:
    input:
        "Enumeration/mapped_reads.dir/mapped/{sample}.bam",
        gtfin
    output:
        "Enumeration/orf_counts.dir/{sample}.tsv.gz"
    threads:
        int(config["Enumerate"]["featureCounts"]["threads"])
    resources:
        mem_mb=int(config["Enumerate"]["featureCounts"]["memory"])
    run:
        command=F.FeatureCounts(input[0],input[1],output[0],config).build()
        shell(command)

#get the features from the config
single_features=config["Enumerate"]["General"]["feature_list"].split(",")
paired_features=[x.replace(":","_BY_")for x in config["Enumerate"]["General"]["feature_pairs"].split(",")]
all_features= single_features+paired_features

#count the chosen annotations by summing the ORF counts
rule countannotations:
    input:
        "Enumeration/orf_counts.dir/{sample}.tsv.gz",
        gtfin
    output:
        expand("Enumeration/annotation_counts.dir/{annotation}/{{sample}}.tsv.gz",annotation=single_features)
    threads:
        int(config["Enumerate"]["countAnnotations"]["threads"])
    resources:
        mem_mb=int(config["Enumerate"]["countAnnotations"]["memory"])
    run:
        command="python {}/workflow/scripts/countFeat.py --orf_counts {} --features {} --multimethod {} --gtf {} --outdir Enumeration/annotation_counts.dir/".format(workflow.basedir,input[0],config["Enumerate"]["General"]["feature_list"],config["Enumerate"]["General"]["multimethod"],input[1].replace(".short",""))
        shell(command)

#count the features as pairs
rule countpairs:
    input:
        "Enumeration/orf_counts.dir/{sample}.tsv.gz",
        gtfin
    output:
        expand("Enumeration/annotation_counts.dir/{annotation}/{{sample}}.tsv.gz",annotation=paired_features)
    threads:
        int(config["Enumerate"]["countAnnotations"]["threads"])
    resources:
        mem_mb=int(config["Enumerate"]["countAnnotations"]["memory"])
    run:
        command="python {}/workflow/scripts/countFeatPairs.py --orf_counts {} --feature_pairs {} --multimethod {} --gtf {} --outdir Enumeration/annotation_counts.dir/".format(workflow.basedir,input[0],config["Enumerate"]["General"]["feature_pairs"],config["Enumerate"]["General"]["multimethod"],input[1].replace(".short",""))
        shell(command)

#combine the counts across samples
rule combinecounts:
    input:
        expand("Enumeration/annotation_counts.dir/{{annotation}}/{sample}.tsv.gz",sample=samples)
    output:
        "Enumeration/combined_counts.dir/{annotation}.tsv"
    threads:
        int(config["Enumerate"]["combineCounts"]["threads"])
    resources:
        mem_mb=int(config["Enumerate"]["combineCounts"]["memory"])
    shell:
        "python {workflow.basedir}/workflow/scripts/combineCounts.py --feature {wildcards.annotation} --countdir Enumeration/annotation_counts.dir/{wildcards.annotation} --outfile {output}"

#generate the enumeration report
rule enumeratereport:
    input:
        expand("Enumeration/combined_counts.dir/{annotation}.tsv",annotation=all_features)
    output:
        report("Enumeration/enumerate_report.html", category="Enumeration")
    params:
        features=",".join(single_features),
        cwd=cwd
    threads:
        int(config["Enumerate"]["report"]["threads"])
    resources:
        mem_mb=int(config["Enumerate"]["report"]["memory"])
    shell:
        'R -e "rmarkdown::render(\'{workflow.basedir}/workflow/report/enumeration_report.Rmd\',output_file=\'{params.cwd}/{output}\')" --args {params.features} {params.cwd}/Enumeration/combined_counts.dir'
