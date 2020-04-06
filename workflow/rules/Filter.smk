#Rules to filter both/one of host reads and rRNA reads from input data
#function to generate the call to sortmerna
from workflow.functions import SortMeRNA as S
from workflow.functions import Hisat2 as H
from workflow.functions import ShellFunctions as SF
from workflow.functions import CountFiltered as CF

configfile: workflow.basedir+"/config/filter.yaml"
            
#run rrna filtering
rule runsortmerna:
    input:
        "{sample}.{ext}"
    output:
        "Filtering/rrna_filter_out.dir/non_rrna/{sample}.{ext}"
    threads:
        int(config["Filter"]["SortMeRNA"]["threads"])
    resources:
        mem=int(config["Filter"]["SortMeRNA"]["memory"])*1000
    run:
        if config["Filter"]["General"]["rrna_filter"] == "true":
            #align to rrnas using sortmerna
            command=S.SortMeRNA(input[0],output[0],config["Filter"]).build()
        else:
            command=SF.symlink(input[0],output[0])
        shell(command)
        
#run genome filtering
#ensure unmapped reads are kept
config["Filter"]["HISAT2"]["keep_unmapped"] == "true"
rule runhisatfilter:
    input:
        "Filtering/rrna_filter_out.dir/non_rrna/{sample}.{ext}"
    output:
        "Filtering/genome_filter_out.dir/unmapped/{sample}.{ext}"
    threads:
        int(config["Filter"]["HISAT2"]["threads"])
    run:
        if config["Filter"]["General"]["host_filter"] == "true":
            #map using Hisat2
            command=H.Hisat2(input[0],output[0],config["Filter"]).build()
        else:
            command=SF.symlink(input[0],output[0])
        shell(command)    


#summarise the counts in each file
rule summarisecounts:
    input:
        orig="{sample}.{ext}",
        rrna="Filtering/rrna_filter_out.dir/non_rrna/{sample}.{ext}",
        host="Filtering/genome_filter_out.dir/unmapped/{sample}.{ext}"
    output:
        "Filtering/report.dir/per_file_summaries/{sample}.{ext}"
    run:
        command=CF.CountFiltered(input[0],input[1],input[2],output,config)
        shell(command)
        
#generate the filtering report
rule filterreport:
    input:
        expand("Filtering/report.dir/per_file_summaries/{sample}.{ext}",zip, sample=samples, ext=exts)
    output:
        "Filtering/report.dir/combined.filter.summary.txt",
        report("Filtering/filter_report.html", category="Filtering")
    run:
        combinedcounts = open(output[0],"w")
        combinedcounts.write("File\tInput\tPost_rRNA_Filtering\tPost_Genome_Filtering\n")
        for i in input:
            sumfile = open(i,"rU").readlines()
            combinedcounts.write(sumfile[1])
        combinedcounts.close()
        command = 'R -e "rmarkdown::render(\'{workflow.basedir}/workflow/report/filter_report.Rmd\',output_file=\'{cwd}/Filtering/filter_report.html\')" --args {cwd}/Filtering/report.dir/combined.filter.summary.txt'
        shell(command)
        
