#Rules to assemble contigs from short reads using metagenomic focused methods
from workflow.functions import Assembler as A
from workflow.functions import ShellFunctions as SF

configfile: workflow.basedir+"/config/assemble.yaml"

#get the tool
Assembler=config["Assemble"]["General"]["assembler"]

#function to allow switching between assembling only or following filtering
def assemblyin(wc):
    prefix="Filtering/genome_filter_out.dir/unmapped/"
    if config["Global"]["run"]=="assemble":
        prefix="./"
    return(prefix+"{{sample}}.{}".format(fileexts[wc.sample]))
    
    
#run rrna filteringxs
rule runassembly:
    input:
        assemblyin
    output:
        "Assembly/contigs.dir/{sample}.contigs.fa.gz"
    threads:
        int(config["Assemble"][Assembler]["threads"])
    resources:
        mem=int(config["Assemble"][Assembler]["memory"])*1000
    run:
        shell("mkdir -p Assembly/assembler.dir/")
        command=A.Assembler(input[0],output[0],Assembler,config["Assemble"][Assembler]).build()
        shell(command)
        
#summarise contigs
rule summarisecontigs:
    input:
        "Assembly/contigs.dir/{sample}.contigs.fa.gz"
    output:
        "Assembly/report.dir/contig_summaries/{sample}.contigs.summary"
    shell:
        "stats.sh in={input} extended=t format=5 shist={output}.shist gchist={output}.gchist overwrite=true > {output}"

#merge summaries and generate the report
rule assemblyreport:
    input:
        expand("Assembly/report.dir/contig_summaries/{sample}.contigs.summary",sample=samples)
    output:
        "Assembly/report.dir/combined.contigs.summary",
        report("Assembly/assemble_report.html", category="Assembly")
    run:
        command=SF.mergecontigsumm(input,output[0],config)
        shell(command)
        command = 'R -e "rmarkdown::render(\'{workflow.basedir}/workflow/report/assemble_report.Rmd\',output_file=\'{cwd}/Assembly/assemble_report.html\')" --args {cwd}/Assembly/report.dir/combined.contigs.summary'
        shell(command)





        
