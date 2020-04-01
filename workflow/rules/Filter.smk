#Rules to filter both/one of host reads and rRNA reads from input data
#function to generate the call to sortmerna
from workflow.functions import SortMeRNA as S
from workflow.functions import Hisat2 as H
from workflow.functions import ShellFunctions as SF

configfile: workflow.basedir+"/config/filter.yaml"
            
#run rrna filtering
rule runsortmerna:
    input:
        "{file}"
    output:
        "Filtering/rrna_filter_out.dir/non_rrna/{file}"
    threads:
        int(config["Filter"]["SortMeRNA"]["threads"])
    resources:
        mem=int(config["Filter"]["SortMeRNA"]["memory"])*1000
    run:
        for i in range(len(input)):
            if config["Filter"]["General"]["rrna_filter"] == "true":
                #align to rrnas using sortmerna
                command=S.SortMeRNA(input[i],output[i],config["Filter"]).build()
            else:
                command=SF.symlink(input[i],output[i])
            shell(command)
        
#run genome filtering
rule runhisatfilter:
    input:
        "Filtering/rrna_filter_out.dir/non_rrna/{file}"
    output:
        "Filtering/genome_filter_out.dir/unmapped/{file}"
    threads:
        int(config["Filter"]["HISAT2"]["threads"])
    run:
        #ensure unmapped reads are kept
        config["Filter"]["HISAT2"]["keep_unmapped"] == "true"
        for i in range(len(input)):
            if config["Filter"]["General"]["host_filter"] == "true":
                #map using Hisat2
                command=H.Hisat2(input[i],output[i],config["Filter"]).build()
            else:
                command=SF.symlink(input[i],output[i])
            shell(command)    
