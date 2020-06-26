'''
object for sequencing data
    
On initialisation:
checks the format from the name
checks for compression from name
check for paired-endness 
'''

import os
import re
import glob
import gzip

#function to open files that checks for gz
def open_file(filename):
    if re.search("gz$",filename):
        return gzip.open(filename,"rt",encoding="utf-8")
    else:
        return open(filename,"r",encoding="utf-8")

class SequencingData:
    def __init__(self,infile):
        self.seqfile_regex =r"(\S+).(fasta$|fasta.gz$|1.fasta.gz$|1.fasta$|fna$|fna.gz$|1.fna.gz$|1.fna$|fa$|fa.gz$|1.fa.gz$|1.fa$|fastq$|fastq.gz$|1.fastq.gz$|1.fastq$|fq$|fq.gz$|1.fq.gz$|1.fq$|contigs.fa.gz$)"
        self.r2_regex =r"(\S+).(2.fasta.gz$|2.fasta$|2.fna.gz$|2.fna$|2.fa.gz$|2.fa$|2.fastq.gz$|2.fastq$|2.fq.gz$|2.fq$)"
        self.r1_regex =r"(\S+).(1.fasta.gz$|1.fasta$|1.fna.gz$|1.fna$|1.fa.gz$|1.fa$|1.fastq.gz$|1.fastq$|1.fq.gz$|1.fq$)"
        self.fileformat =  None
        self.paired = False
        self.interleaved = False
        self.compressed = False
        self.pairedname = ""
        self.filepath = infile
        self.filename = os.path.basename(infile)
        self.openfile = None
        self.cleanname = None
        self.head = []
        #check file can be opened on init & capture header (first 5 lines used for interleave and format checks)
        try:
            self.openfile = open_file(self.filepath)
        except FileNotFoundError as e:
            raise Exception("cannot open file {}".format(self.filepath)) from e
        self.head = [self.openfile.readline().rstrip("\n") for x in range(5)]
        #autocheck file format and pairedness, read count must be specified seperately
        self.getFormat()
        self.isPaired()
        self.cleanName()
        
    #check it is fasta or fastq and if compressed    
    def getFormat(self):
        assert re.search(self.seqfile_regex,self.filepath),"file {} is not of the correct format (fasta or fastq) or has incorrect file name structure.".format(self.filename)
        assert not(re.search(self.r2_regex,self.filepath)), "Read 2 file provided ({}) please use read 1 file".format(self.filename)
        if re.search("(.fq|.fastq)",self.filepath):
            self.fileformat="fastq"
        else:
            self.fileformat="fasta"
        if self.filepath.endswith(".gz"):
            self.compressed = True
        if self.fileformat == "fasta":
            assert self.head[0][0] == ">", "invalid header on first line for fasta format"
        else:
            assert self.head[0][0] == "@", "invalid header on first line for fastq format"

            
    #check if paired and if containts interleaved pairs or matching files
    def isPaired(self):
        if re.search(self.r1_regex,self.filepath):
            paired_name = self.filepath.replace(".1.",".2.")
            assert len(glob.glob(paired_name)) > 0, "cannot find read 2 file at location {} associated with read 1 file {}".format(paired_name,self.filename)
            paired_name = os.path.basename(paired_name)
            self.paired = True
            self.pairedname = paired_name
        elif self.isInterleaved(self.filepath):
            self.paired = True
            self.interleaved = True

            
    #abstract out the interleaving check for clarity
    def isInterleaved(self,fpath):
        pairindex = 2
        if self.fileformat == "fastq":
            pairindex = 4
        if self.head[0].endswith("/1"):
            assert self.head[0].strip("/1") == self.head[pairindex].strip("/2"), "first read in file {} is named in the interleaved format (/1) but does not have a matching read 2 as expected".format(self.filename)
            return True
        else:
            return False

    #name with no extension to use later
    def cleanName(self):
        if self.paired == True:
            self.cleanname=re.search(self.r1_regex,self.filename).group(1)
        else:
            self.cleanname=re.search(self.seqfile_regex,self.filename).group(1)
