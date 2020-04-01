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
        extensions=("fna","fa","fasta","fastq")
        for i in extensions:
            assert not(self.filepath.endswith((".2",".2.gz"))), "Read 2 file provided ({}) please use read 1 file".format(self.filename)
            if self.filepath.endswith((i,i+".1.gz",i+".gz",i+".1")):
                if i == "fastq":
                    self.fileformat=i
                else:
                    self.fileformat="fasta"
        if self.filepath.find(".gz")!=-1:
            self.compressed = True
        assert self.fileformat,"file {} is not of the correct format (fasta or fastq).".format(self.filename)
        if self.fileformat == "fasta":
            assert self.head[0][0] == ">", "invalid header on first line for fasta format"
        else:
            assert self.head[0][0] == "@", "invalid header on first line for fastq format"

            
    #check if paired and if containts interleaved pairs or matching files
    def isPaired(self):
        if self.filepath.endswith((".1",".1.gz")):
            paired_name = self.filepath.replace(".1",".2")
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
        seqfile_regex =r"(\S+).(fasta$|fasta.gz|fasta.1.gz|fasta.1|fna$|fna.gz|fna.1.gz|fna.1|fa$|fa.gz|fa.1.gz|fa.1|fastq$|fastq.gz|fastq.1.gz|fastq.1)"
        self.cleanname=re.search(seqfile_regex,self.filename).group(1)
