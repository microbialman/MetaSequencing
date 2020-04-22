from argparse import ArgumentParser
import re, glob, gzip

#get the input and output files
parser = ArgumentParser()
parser.add_argument("--feature", dest="feat", help="feature to combine counts for")
parser.add_argument("--countdir", dest="cdir", help="head directory contianing sub-directories with feature counts")
parser.add_argument("--outfile",dest="outfile", help="output destination for combined table")

args = parser.parse_args()

#get all the files in the count directory
sampfiles=glob.glob("{}/*.tsv.gz".format(args.cdir))

#dictionary to hold the counts
featcounts={}
samplenames=[]
featnames=[]

#parse the files
for i in range(len(sampfiles)):
    fileopen=gzip.open(sampfiles[i],'rt').readlines()
    sample=re.search("{}/(\S+).tsv.gz".format(args.cdir),sampfiles[i]).group(1)
    samplenames.append(sample)
    for j in fileopen[1:]:
        row=j.strip("\n").split("\t")
        if row[0] not in featcounts:
            featnames.append(row[0])
            featcounts[row[0]]=[0.0]*len(sampfiles)
        featcounts[row[0]][i]=float(row[1])

#write the output
outfile=open(args.outfile,"w")
outfile.write("{}\t{}\n".format(args.feat,"\t".join(samplenames)))
for i in featnames:
    outfile.write("{}\t{}\n".format(i,"\t".join([str(x) for x in featcounts[i]])))
outfile.close()
