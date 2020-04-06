#Rules for functional and taxonomic annotation of contigs from metagenomic data
from workflow.functions import Prodigal as P

configfile: workflow.basedir+"/config/annotate.yaml"

#function to allow switching between assembling only or following filtering
def annotatein(wc):
    prefix="Assembly/contigs.dir/"
    if config["Global"]["run"]=="annotate":
        return("./{{sample}}.{}".format(fileexts[wc.sample]))
    else:
        return("Assembly/contigs.dir/{sample}.contigs.fa.gz")
    
    
#detect ORFs using prodigal
rule runprodigal:
    input:
        annotatein
    output:
        "Annotation/orfs.dir/{sample}.orf_peptides.gz"
    threads:
        int(config["Annotate"]["Prodigal"]["threads"])
    resources:
        mem=int(config["Annotate"]["Prodigal"]["memory"])*1000
    run:
        pobj=P.Prodigal(input[0],output[0],config)
        shell(pobj.build())

#split the orf files to chunks to increase parallelisation
rule chunkorfs:
    input:
        "Annotation/orfs.dir/{sample}.orf_peptides.gz"
    output:
        "Annotation/functional_annotations.dir/emapper_chunks/{sample}.chunk.log"
    run:
        shell("python {}/workflow/scripts/fastaToChunks.py --input {} --output_prefix {} --chunk_size {}".format(
            workflow.basedir,input[0],output[0].replace(".chunk.log",""),config["Annotate"]["Eggnogmapper"]["chunksize"]))
