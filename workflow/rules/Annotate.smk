#Rules for functional and taxonomic annotation of contigs from metagenomic data
from workflow.functions import Prodigal as P
from workflow.functions import Eggnogmapper as E
from workflow.functions import Kraken2 as K
import re,os

configfile: workflow.basedir+"/config/annotate.yaml"

#function to allow switching between annotation only or following assembly in whole run
def annotatein(wc):
    prefix="Assembly/contigs.dir/"
    if config["Global"]["run"]=="annotate":
        return("{{sample}}.{}".format(fileexts[wc.sample]))
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
        mem_mb=int(config["Annotate"]["Prodigal"]["memory"])
    run:
        pobj=P.Prodigal(input[0],output[0],config)
        shell(pobj.build())

#split the orf files to chunks to increase parallelisation
checkpoint chunkorfs:
    input:
        "Annotation/orfs.dir/{sample}.orf_peptides.gz"
    output:
        directory("Annotation/functional_annotations.dir/emapper_chunks/{sample}")
    params:
        nreads=config["Annotate"]["Eggnogmapper"]["chunksize"]
    shell:
        "mkdir -p {output} && python {workflow.basedir}/workflow/scripts/fastaToChunks.py --input {input} --outdir {output} --nreads {params.nreads}"
        
#find seed orthologs for each chunk
rule funseed:
    input:
        "Annotation/functional_annotations.dir/emapper_chunks/{sample}/{i}.fa"
    output:
        temp("Annotation/functional_annotations.dir/seed_orthologs/{sample}/{i}.emapper.seed_orthologs")
    threads:
        int(config["Annotate"]["Eggnogmapper"]["threads_orth"])
    resources:
        mem_mb=int(config["Annotate"]["Eggnogmapper"]["memory_orth"])
    run:
        eobj=E.Eggnogmapper(input[0],output[0],config,"ortho")
        shell(eobj.build())
        
#annotate seeds
rule annotseed:
    input:
        "Annotation/functional_annotations.dir/seed_orthologs/{sample}/{i}.emapper.seed_orthologs"
    output:
        temp("Annotation/functional_annotations.dir/chunk_annots/{sample}/{i}.emapper.annotations")
    threads:
        int(config["Annotate"]["Eggnogmapper"]["threads_annot"])
    resources:
        mem_mb=int(config["Annotate"]["Eggnogmapper"]["memory_annot"])
    run:
        eobj=E.Eggnogmapper(input[0],output[0],config,"annot")
        shell(eobj.build())
        
#function to generate input from split emapper chunks
def aggregatechunksin(wc):
    chkpnt = checkpoints.chunkorfs.get(sample=wc.sample).output[0]
    return(expand("Annotation/functional_annotations.dir/chunk_annots/{sample}/{i}.emapper.annotations",
                  sample=wc.sample,
                  i=glob_wildcards(os.path.join(chkpnt,"{i}.fa")).i))
       
#aggregate split functional annotation files
rule aggregatechunks:
    input:
        aggregatechunksin
    output:
        "Annotation/functional_annotations.dir/{sample}.functional.annotations.gz"
    run:
        sample = re.search("functional_annotations\.dir/chunk_annots/(\S+)/[0-9]*\.emapper\.annotations",input[0]).group(1)
        chunkfolder = "Annotation/functional_annotations.dir/emapper_chunks/"+sample
        command = "cat {0} > {1} && gzip {1} && rm -r {2}".format(" ".join(input),output[0].strip(".gz"),chunkfolder)
        shell(command)


#use kraken to assign taxonomy to the contigs
rule runkraken:
    input:
        annotatein
    output:
        temp("Annotation/taxonomic_annotations.dir/kraken_out/{sample}.kraken_out")
    threads:
        int(config["Annotate"]["Kraken2"]["threads"])
    resources:
        mem_mb=int(config["Annotate"]["Kraken2"]["memory"])
    run:
        kobj = K.Kraken2(input[0],output[0],config)
        shell(kobj.callkraken())

#translate kraken ids output to full taxonomic names using TaxonKit
rule translatekraken:
    input:
        "Annotation/taxonomic_annotations.dir/kraken_out/{sample}.kraken_out"
    params:
        datadir=config["Annotate"]["TaxonKit"]["datadir"],
        taxonpre=config["Annotate"]["TaxonKit"]["taxonpre"]
    output:
        "Annotation/taxonomic_annotations.dir/{sample}.taxonomic.annotations.gz"
    shell:
        "python {workflow.basedir}/workflow/scripts/translateKraken2.py --krakenout {input} --translatedout {output} --taxdatadir {params.datadir} --taxaprefixes {params.taxonpre}"

#combine the taxonomic and functional annotations to gtf file
rule combineannotations:
    input:
        orfs="Annotation/orfs.dir/{sample}.orf_peptides.gz",
        fun="Annotation/functional_annotations.dir/{sample}.functional.annotations.gz",
        tax="Annotation/taxonomic_annotations.dir/{sample}.taxonomic.annotations.gz"
    output:
        full="Annotation/combined_annotations.dir/{sample}.orf_annotations.gtf.gz",
        short="Annotation/combined_annotations.dir/{sample}.orf_annotations.short.gtf.gz"
    resources:
        mem_mb=int(config["Annotate"]["Merge"]["memory"])
    shell:
        "python {workflow.basedir}/workflow/scripts/makeGtf.py --orfs {input.orfs} --functions {input.fun} --taxonomy {input.tax} --output {output.full}"
        
#summarise the gtfs for the report
rule summarisegtfs:
    input:
        expand("Annotation/combined_annotations.dir/{sample}.orf_annotations.gtf.gz",sample=samples)
    output:
        "Annotation/report.dir/annotation_summary.tsv.gz",
        "Annotation/report.dir/orf_summary.tsv.gz"
    resources:
        mem_mb=int(config["Annotate"]["Merge"]["memory"])
    run:
        command="python {}/workflow/scripts/annotationSummaryTable.py --gtfs {} --annot-out {} --orf-output {}".format(workflow.basedir,",".join(input),output[0],output[1])
        shell(command)

#get the parings to generate the gmt files for
gmtpairs = [x.replace(":","_TO_") for x in config["Annotate"]["Gmt"]["pairs"].split(",")]
        
#generate gmt files for use in set ernichment analyses
rule makegmts:
    input:
        expand("Annotation/combined_annotations.dir/{sample}.orf_annotations.gtf.gz",sample=samples)
    output:
        "Annotation/gmt_files.dir/{pairing}.gmt.gz"
    resources:
        mem_mb=int(config["Annotate"]["Gmt"]["memory"])
    run:
        infiles=",".join(input)
        pair=wildcards.pairing.split("_TO_")
        command="python {}/workflow/scripts/makeGmt.py --set {} --identifier {} --gtfs {} --outfile {}".format(workflow.basedir,pair[0],pair[1],infiles,output)
        shell(command)

#generate the taxon gmts
rule maketaxongmts:
    input:
        expand("Annotation/combined_annotations.dir/{sample}.orf_annotations.gtf.gz",sample=samples)
    output:
        "Annotation/taxon_gmt_files.dir/species.gmt.gz"
    resources:
        mem_mb=int(config["Annotate"]["Gmt"]["memory"])
    run:
        infiles=",".join(input)
        command="python {}/workflow/scripts/makeTaxonGmt.py --gtfs {} --outdir {}".format(workflow.basedir,infiles,os.path.dirname(str(output)))
        shell(command)
                
#make the annotation report
rule annotationreport:
    input:
        "Annotation/report.dir/annotation_summary.tsv.gz",
        "Annotation/report.dir/orf_summary.tsv.gz",
        expand("Annotation/gmt_files.dir/{pairing}.gmt.gz",pairing=gmtpairs),
        "Annotation/taxon_gmt_files.dir/species.gmt.gz"
    output:
        report("Annotation/annotate_report.html", category="Annotate")
    resources:
        mem_mb=int(config["Annotate"]["Report"]["memory"])
    shell:
        'R -e "rmarkdown::render(\'{workflow.basedir}/workflow/report/annotate_report.Rmd\',output_file=\'{cwd}/Annotation/annotate_report.html\')" --args {cwd}/Annotation/report.dir/annotation_summary.tsv.gz {cwd}/Annotation/report.dir/orf_summary.tsv.gz'
