
#this scrpt generates a summary table of the annotations that can then be used in the report scrit
from argparse import ArgumentParser
import re, glob
import gzip
import numpy

#get the orfs, function and taxonomy annotations from commadn line
parser = ArgumentParser()
parser.add_argument("--gtfs", dest="gtfs", help="List of GTFs containing the combined functional and taxonomic annotations for each sample")
parser.add_argument("--annot-output", dest="outfile", help="Annotation summary output filename")
parser.add_argument("--orf-output", dest="orfout", help="ORF summary output filename")
args = parser.parse_args()

#get the list of gtfs (full not short) for each sample
gtfs = args.gtfs.split(",")

featdic={}
orfdic={}

taxlevs=["kingdom","phylum","class","order","family","genus","species"]
nonlistfun=["Predicted_protein_name"]
#list funs will be split on commas
listfun=["Gene_Ontology_terms","EC_number","KEGG_ko","KEGG_Pathway","KEGG_Module","KEGG_Reaction","BRITE","BiGG_Reaction","eggNOG_OGs","COG_Functional_Category"]

#open the orf output file
orffile=gzip.open(args.orfout,"wt")
orffile.write("Sample\tNo_ORFs\tMean_ORF_Size\tSD_ORF_Size\tMean_ORFS_Per_Contig\tSD_ORFS_Per_Contig\n")

#function to parse features, some are listed (i.e. one annotation conatins many feautes, eg. GO terms)
def featAdd(sample,feat,listed):
    def addDic(key,val):
        if key not in featdic:
            featdic[key]={"type":val,"samples":[0]*len(gtfs)}
        if featdic[key]["samples"][i]==0:
            featdic[key]["samples"][i]=1
    if listed == False:
        if feat[1] != "":
            addDic(feat[1],feat[0])
    else:
        for x in feat[1].split(","):
            if x != "":
                addDic(x,feat[0])

#go through gtf files and record occurences of each feature
totalsizes=[]
totalorfper=[]
totalcount=0
for i in range(len(gtfs)):
    #open annotation
    gtf=gzip.open(gtfs[i],"rt")
    contig_ORF_counts={}
    ORF_sizes=[]
    no_ORFs=0
    #parse
    for j in gtf:
        row=j.strip("\n").split("\t")
        #parse the orf
        contig=row[0]
        length=int(row[4])-int(row[3])
        ORF_sizes.append(length)
        no_ORFs+=1
        if contig in contig_ORF_counts:
            contig_ORF_counts[contig]+=1
        else:
            contig_ORF_counts[contig]=1
        #parse the annotations
        annotations=[x.split('"') for x in row[-1].split(";")]
        for k in annotations[1:]:
            k[0]=k[0].strip(" ")
            if k[0] in nonlistfun or k[0] in taxlevs:
                featAdd(i,k,False)
            elif k[0] in listfun:
                featAdd(i,k,True)
    #write out ORF summary for sample
    samp=gtfs[i]
    noorfs=str(no_ORFs)
    meansize=str(numpy.mean(ORF_sizes))
    sdsize=str(numpy.std(ORF_sizes))
    orfspercon=[int(x) for x in contig_ORF_counts.values()]
    meanorfspercon=str(numpy.mean(orfspercon))
    sdorfspercon=str(numpy.std(orfspercon))
    orffile.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(samp,noorfs,meansize,sdsize,meanorfspercon,sdorfspercon))
    #update totals
    totalcount+=no_ORFs
    totalsizes+=ORF_sizes
    totalorfper+=orfspercon
    gtf.close()
#write totals
orffile.write("TOTAL\t{}\t{}\t{}\t{}\t{}\n".format(str(totalcount),str(numpy.mean(totalsizes)),str(numpy.std(totalsizes)),str(numpy.mean(totalorfper)),str(numpy.std(totalorfper))))
orffile.close()


#write the annotation output file
outfile=gzip.open(args.outfile,"wt")
outfile.write("feat_name\tfeat_type\t"+"\t".join(gtfs)+"\n")
for i in featdic:
    outfile.write("{}\t{}\t{}\n".format(i,featdic[i]["type"],"\t".join([str(x) for x in featdic[i]["samples"]])))
outfile.close()
