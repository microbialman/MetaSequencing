from workflow.functions import SequencingData

'''
class to build call to featureCounts
'''

class FeatureCounts:
    def __init__(self,samfile,gtffile,outfile,params):
        self.sam = samfile
        self.gtf = gtffile
        self.out = outfile.replace(".gz","")
        self.params = params
        self.statementlist = []
        self.buildStatement()
        self.zipOut()
        
    def buildStatement(self):
        flist=["featureCounts"]
        flist.append("-a {} -o {} -t ORF -f".format(self.gtf,self.out))
        if self.params["Enumerate"]["featureCounts"]["additional_args"] != "":
            flist.append(self.params["Enumerate"]["featureCounts"]["additional_args"])
        flist.append(self.sam)
        self.statementlist.append(" ".join(flist))

    def zipOut(self):
        self.statementlist.append("gzip {}".format(self.out))

    def build(self):
        return(" && ".join(self.statementlist))
