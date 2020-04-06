import os
from workflow.functions import SequencingData


'''
general class for calls to assembly tools
'''

class Assembler():
    def __init__(self,seqdat,outfinal,tool,params):
        self.seqdat = SequencingData.SequencingData(seqdat)
        self.f1 = self.seqdat.filepath
        if self.seqdat.paired == True:
            self.f2 = os.path.dirname(self.seqdat.filepath)+"/"+self.seqdat.pairedname
        self.outfinal = outfinal
        self.tool = tool
        self.params = params
        self.outfile = ""
        self.statementlist = []
        self.cdir = os.getcwd()
        if self.tool == "MEGAHIT":
            self.megahit()
        elif self.tool == "metaSPAdes":
            self.metaspades()
        self.linker()
        
    #functions to build calls for each assembler    
    def megahit(self):
        mcall = ["megahit"]
        if self.seqdat.paired == True:
            if self.seqdat.interleaved == True:
                mcall.append("--12 {}/{}".format(self.cdir,self.f1))
            else:
                mcall.append("-1 {0}/{1} -2 {0}/{2}".format(self.cdir,self.f1,self.f2))
        else:
            mcall.append("-r {}/{}".format(self.cdir,self.f1))
        if self.params["threads"] != "":
            mcall.append("-t {}".format(self.params["threads"]))
        mcall.append("-o {}/Assembly/assembler.dir/{}".format(self.cdir,self.seqdat.cleanname))
        mcall.append("--out-prefix {}".format(self.seqdat.cleanname))
        if self.params["additional_args"] != "":
            mcall.append(self.params["additional_args"])
        self.statementlist.append(" ".join(mcall))
        #set output and remove intermediate contigs
        self.outfile="{0}/Assembly/assembler.dir/{1}/{1}.contigs.fa".format(self.cdir,self.seqdat.cleanname)
        self.statementlist.append("rm -r {}/Assembly/assembler.dir/{}/intermediate_contigs".format(self.cdir,self.seqdat.cleanname))
    
    def metaspades(self):
        mcall = ["metaspades.py"]
        if self.seqdat.paired == True:
            if self.seqdat.interleaved == True:
                mcall.append("--12 {}/{}".format(self.cdir,self.f1))
            else:
                mcall.append("-1 {0}/{1} -2 {0}/{2}".format(self.cdir,self.f1,self.f2))
        else:
            mcall.append("-s {}/{}".format(self.cdir,self.f1))
        mcall.append("-o {}/Assembly/assembler.dir/{}".format(self.cdir,self.seqdat.cleanname))
        if self.params["threads"] != "":
            mcall.append("-t {}".format(self.params["threads"]))
        if self.params["memory"] != "":
            mcall.append("-m {}".format(self.params["memory"]))
        if self.params["additional_args"] != "":
            mcall.append(self.params["additional_args"])
        self.statementlist.append(" ".join(mcall))
        #set output and remove intermediate contigs
        self.outfile="{0}/Assembly/assembler.dir/{1}/contigs.fasta".format(self.cdir,self.seqdat.cleanname)
        #self.statementlist.append("rm -r {}/Assembly/assembler/{}/intermediate_contigs".format(self.cdir,self.seqdat.cleanname))

    def linker(self):
        #zip and link to contig folder
        self.statementlist.append("gzip {}".format(self.outfile))
        self.statementlist.append("ln -rs {1}.gz {0}/{2}".format(self.cdir,self.outfile,self.outfinal))
    
    def build(self):
        statement = " && ".join(self.statementlist)
        return(statement)
