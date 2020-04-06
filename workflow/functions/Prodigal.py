from workflow.functions import SequencingData
import os

'''
class to build call to Prodigal
'''
class Prodigal:    
    def __init__(self,seqdat,outfile,params):
        
        self.seqdat = SequencingData.SequencingData(seqdat)
        self.outfile = outfile
        self.params = params
        self.cat="cat"
        self.short=outfile.replace("_peptides.gz","")
        self.statementlist = []
        self.checkcompress()
        self.callprodigal()
        self.compress()
        
    def checkcompress(self):
        if self.seqdat.compressed==True:
            self.cat="zcat"

    def callprodigal(self):
        plist=[]
        plist.append("{} {} | prodigal -o {} -a {} -p meta -q".format(
            self.cat,self.seqdat.filepath,self.short+"_positions",self.short+"_peptides"))
        if self.params["Annotate"]["Prodigal"]["additional_args"] != "":
            plist.append(self.params["Annotate"]["Prodigal"]["additional_args"])
        self.statementlist.append(" ".join(plist))

    def compress(self):
        self.statementlist.append("gzip {0}_positions && gzip {0}_peptides".format(self.short))
        
    def build(self):
        return(" && ".join(self.statementlist))
