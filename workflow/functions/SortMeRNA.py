import os
import sys

from workflow.functions import SequencingData

'''
class to build call to SortMeRNA
'''

class SortMeRNA:    
    def __init__(self,seqdat,outfile,params):
        
        self.seqdat = SequencingData.SequencingData(seqdat)
        self.outfile = outfile
        self.outdir = os.getcwd()+"/"+os.path.dirname(outfile).replace("/non_rrna","")
        self.indir = os.getcwd()+"/"
        self.params = params
        self.filelocation = ""
        self.statementlist = []
        self.rminter = False
        self.checkInterleave()
        self.buildStatement()
        self.deInterleave()
        self.deletemapped()
        
    #add statements to generate interleaved files if necessary
    def checkInterleave(self):
        #if single end or already interleaved just use the original file location
        if self.seqdat.paired == False or self.seqdat.interleaved == True:
            self.filelocation = os.getcwd()+"/"+self.seqdat.filepath
        #else generate commands to make temp interleaved reads
        else:
            self.statementlist.append("mkdir -p {}/interleaved && mkdir -p {}/rrna".format(self.outdir,self.outdir))
            self.filelocation = self.outdir+"/interleaved/{}".format(self.seqdat.cleanname+"."+self.seqdat.fileformat)
            self.statementlist.append("seqtk mergepe {} {} >{}".format(self.indir+self.seqdat.filepath,self.indir+self.seqdat.pairedname,self.filelocation))

    #make the main call to sortmerna
    def buildStatement(self):
        sortlist = ["sortmerna"]
        sortlist.append(self.refList())
        if self.params["SortMeRNA"]["paired"] == "out":
            sortlist.append("--paired_out")
        else:
            sortlist.append("--paired_in")
        sortlist.append("--reads {}".format(self.filelocation))
        sortlist.append("--aligned {}".format(self.outdir+"/rrna/"+self.seqdat.cleanname))
        sortlist.append("--other {}".format(self.outdir+"/non_rrna/"+self.seqdat.cleanname))
        sortlist.append("--fastx")
        if self.params["SortMeRNA"]["additional_args"] != "":
            sortlist.append(self.params["SortMeRNA"]["additional_args"])
        sortlist.append("-a {}".format(self.params["SortMeRNA"]["threads"]))
        if self.params["SortMeRNA"]["memory"] != "false":
            sortlist.append("-m {}".format(str(int(self.params["SortMeRNA"]["memory"])*900*int(self.params["SortMeRNA"]["threads"]))))
        self.statementlist.append(" ".join(sortlist))
        #compressed aligned output
        self.statementlist.append("gzip {}.*".format(self.outdir+"/rrna/"+self.seqdat.cleanname))

    #if input was interleaved remove interleaved temp file (if necessary)
    def deInterleave(self):
        if self.seqdat.paired == False:
            self.statementlist.append("gzip {}.*".format(self.outdir+"/non_rrna/"+self.seqdat.cleanname))
        elif self.seqdat.paired == True:
            if self.seqdat.interleaved == False:
                #remove temp interleave file if one was made
                self.statementlist.append("rm {}".format(self.filelocation))
                #undo interleaving
                otherf = self.outdir+"/non_rrna/"+self.seqdat.cleanname+"."+self.seqdat.fileformat
                pair1 = self.outdir+"/non_rrna/"+self.seqdat.filename
                pair2 = self.outdir+"/non_rrna/"+self.seqdat.pairedname
                self.statementlist.append("seqtk seq -l0 -1 {} > {}".format(otherf,pair1.strip(".gz")))
                self.statementlist.append("seqtk seq -l0 -2 {} > {}".format(otherf,pair2.strip(".gz")))
                self.statementlist.append("rm {}".format(otherf))
                #compress the outputs
                self.statementlist.append("gzip {}".format(pair1.strip(".gz")))
                self.statementlist.append("gzip {}".format(pair2.strip(".gz")))
            else:
                self.statementlist.append("gzip {}.*".format(self.outdir+"/non_rrna/"+self.seqdat.cleanname))
                    
    #abstract out making reference command as it is long 
    def refList(self):
        reffastas = self.params["SortMeRNA"]["ref_fastas"].split(",")
        if self.params["SortMeRNA"]["ref_index"] != "false":
            refindex = self.params["SortMeRNA"]["ref_index"].split(",")
        else:
            refindex = reffastas
            refindex = [os.path.basename(i) for i in refindex]
            refindex = [os.getcwd()+"/ref_index.dir/{}-db".format(i) for i in refindex]
        paired = [",".join([x[0],x[1]]) for x in zip(reffastas,refindex)]
        return("--ref {}".format(":".join(paired)))
            
    def deletemapped(self):
        if self.params["SortMeRNA"]["delete_mapped"] == "true":
            self.statementlist.append("rm {}.*".format(self.outdir+"/rrna/"+self.seqdat.cleanname))
    
    def build(self):
        return(" && ".join(self.statementlist))
