#this script generates GMT files across gtf files for later use in set enrichment approaches
#this version generates taxon set lists at different base resolutions (e.g. all higher level groupings of species ids)
from argparse import ArgumentParser
import gzip

#get the orfs, function and taxonomy annotations from commadn line
parser = ArgumentParser()
parser.add_argument("--gtfs", dest="gtfs", help="Comma seperated list of gtfs to parse gmt from.")
parser.add_argument("--outdir", dest="outdir", help="Output file directory location.")
args = parser.parse_args()

#get the list of gtfs for each sample
gtfs = args.gtfs.split(",")

#dictionary to store the taxa
levels=["kingdom","phylum","class","order","family","genus","species"]
taxdic={x:set() for x in levels}

#get all the taxon info
for i in gtfs:
    if i[-3:] == ".gz":
        fileopen=gzip.open(i,"rt")
    else:
        fileopen=open(i,"r")
    for j in fileopen:
        row=j.strip("\n").strip(";")
        row=row.split("\t")
        annots=[x.split() for x in row[-1].split(";")]
        annotdic={x[0].replace(" ",""):x[1].replace('"','') for x in annots}
        for k in levels:
            if k in annotdic:
                taxdic[k].add(annotdic[k])

#write out the gmt files at each level
for i in range(1,len(levels)):
    vals=taxdic[levels[i]]
    mappingdic={}
    for j in vals:
        split=j.split("|")
        for k in range(len(split)-1):
            tname="|".join(split[0:k+1])
            if tname not in mappingdic:
                mappingdic[tname]=set()
            mappingdic[tname].add(j)
    outfile=gzip.open(args.outdir+"/{}.gmt.gz".format(levels[i]),"wt")
    for j in mappingdic:
        outfile.write("{}\tNA\t{}\n".format(j,"\t".join(mappingdic[j])))
    outfile.close()
