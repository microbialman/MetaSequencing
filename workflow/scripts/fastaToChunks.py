from argparse import ArgumentParser
import re
import gzip

#script to split fastsa file by number of reads (can't use split for multi-line fastas)
parser = ArgumentParser()
parser.add_argument("--input", dest="inputfile", help="input fasta")
parser.add_argument("--outdir", dest="out", help="output directory")
parser.add_argument("--nreads", dest="nreads", help="no. of reads per chunk")
args = parser.parse_args()

#open the input file
if re.search("gz$",args.inputfile):
    fasta=gzip.open(args.inputfile,"rt",encoding="utf-8")
else:
    fasta=open(args.inputfile,"r",encoding="utf-8")

#other vars
odir = args.out
n = int(args.nreads)

#initiate counts
current_chunk=1
read_count=0

#open first chunk file
chunkfile=open(odir+"/{}.fa".format(current_chunk),"w")

#write chunks
for i in fasta:
    if i[0]==">":
        if read_count == n:
            chunkfile.close()
            current_chunk+=1
            chunkfile=open(odir+"/{}.fa".format(current_chunk),"w")
            read_count = 0
        else:
            read_count+=1
    chunkfile.write(i)

chunkfile.close()
