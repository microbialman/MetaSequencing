import re, os
from workflow.functions import SequencingData

'''
class to build call to HISAT2
'''

class Hisat2:    
    def __init__(self,seqdat,outfile,params):
        self.seqdat = SequencingData.SequencingData(seqdat)
        self.outfile = outfile
        self.outclean = "".join(self.outfile.strip(".gz").split(".")[0:-1])
        self.outdir = os.getcwd()+"/"+os.path.dirname(outfile)
        self.outdirmapped = self.outdir.replace("/unmapped","/mapped")
        self.outdirunmapped = self.outdir.replace("/mapped","/unmapped")
        self.indir = os.getcwd()+"/"
        self.params = params
        self.statementlist = []
        self.mkdirs()
        self.buildStatement()
        self.moveUnmapped()
        self.bam()
        self.deletemapped()
        
    def mkdirs(self):
        self.statementlist.append("mkdir -p {}".format(self.outdirmapped))
        if self.params["HISAT2"]["keep_unmapped"] == "true":
            self.statementlist.append("mkdir -p {}".format(self.outdirunmapped))
        
    #make the main call to Hisat
    def buildStatement(self):
        satlist = ["hisat2"]
        satlist.append("-x {}".format(self.params["HISAT2"]["ref_index"]))
        if self.seqdat.fileformat == "fastq":
            satlist.append("-q")
        else:
            satlist.append("-f")
        if self.seqdat.paired == True:
            satlist.append("-1 {0}{1} -2 {0}{2}/{3}".format(self.indir,self.seqdat.filepath,os.path.dirname(self.seqdat.filepath),self.seqdat.pairedname))
        else:
            satlist.append("-U {}".format(self.indir+self.seqdat.filepath))
        satlist.append("--threads {}".format(self.params["HISAT2"]["threads"]))
        satlist.append("-S {}/{}.sam".format(self.outdirmapped,self.seqdat.cleanname))
        if self.params["HISAT2"]["additional_args"] != "":
            satlist.append(self.params["HISAT2"]["additional_args"])
        if self.params["HISAT2"]["keep_unmapped"] == "true":
            if self.seqdat.paired == "false":
                satlist.append("--un-gz {}/{}".format(self.outdirunmapped,self.seqdat.cleanname))
            else:
                satlist.append("--un-conc-gz {}/{}".format(self.outdirunmapped,self.seqdat.cleanname))
        self.statementlist.append(" ".join(satlist))

    #rename the unmapped reads to match original files if keeping them
    def moveUnmapped(self):
        if self.params["HISAT2"]["keep_unmapped"] == "true":
            if self.seqdat.paired == "false":
                self.statementlist.append("mv {0}/{1} {0}/{2}".format(self.outdirunmapped,self.seqdat.cleanname,self.seqdat.filename.replace(".gz","")+".gz"))
            else:
                self.statementlist.append("mv {0}/{1}.1 {0}/{2}".format(self.outdirunmapped,self.seqdat.cleanname,self.seqdat.filename.replace(".gz","")+".gz"))
                self.statementlist.append("mv {0}/{1}.2 {0}/{2}".format(self.outdirunmapped,self.seqdat.cleanname,self.seqdat.pairedname.replace(".gz","")+".gz"))
        
    #convert sam to bam to save space in chosen
    def bam(self):
        if self.params["HISAT2"]["convert_to_bam"] == "true":
            self.statementlist.append("samtools view -bS {0}/{1}.sam > {0}/{1}.bam && rm {0}/{1}.sam".format(self.outdirmapped,self.seqdat.cleanname))

    #delete mapped files if just using for filtering
    def deletemapped(self):
        if self.params["HISAT2"]["delete_mapped"] == "true":
            self.statementlist.append("rm {}/{}.*".format(self.outdirmapped,self.seqdat.cleanname))
            
    def build(self):
        return(" && ".join(self.statementlist))
