#this script generates GMT files across gtf files for later use in set ernichment approaches
from argparse import ArgumentParser
import gzip

#get the orfs, function and taxonomy annotations from commadn line
parser = ArgumentParser()
parser.add_argument("--set", dest="st", help="Identifier defining a set (e.g. KEGG pathway)")
parser.add_argument("--identifier", dest="iden", help="Identfier to be allocated to sets (e.g. protein name).")
parser.add_argument("--gtfs", dest="gtfs", help="Comma seperated list of gtfs to parse gmt from.")
parser.add_argument("--outfile", dest="outfile", help="Output file location.")
args = parser.parse_args()

#get the identifiers from args
st=args.st
iden=args.iden


#get the list of gtfs for each sample
gtfs = args.gtfs.split(",")

#parse all the files and build a dictionary of the sets
setdic={}

#parse the gtfs 
for i in gtfs:
    fileopen=gzip.open(i,"rt")
    for j in fileopen:
        annots=[x.split() for x in j.split("\t")[8].split(";")]
        if len(annots)==2:
            pass
        else:
            rdic=dict(annots)
            if st in rdic and iden in rdic:
                sets=rdic[st].replace('"','')
                idens=rdic[iden].replace('"','')
                if sets=="" or idens=="":
                    pass
                else:
                    for k in sets.split(","):
                        if k not in setdic:
                            setdic[k]=set()
                        for l in idens.split(","):
                            setdic[k].add(l)
    fileopen.close()

#write to the output file
output=gzip.open(args.outfile,"wt")
for i in setdic:
    output.write("{}\tNA\t{}\n".format(i,"\t".join(setdic[i])))
output.close()
