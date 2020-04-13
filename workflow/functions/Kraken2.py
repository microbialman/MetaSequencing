'''
class to build calls to Kraken2 and translating script
'''
class Kraken2:    
    def __init__(self,infile,outfile,params):
        
        self.infile = infile
        self.outfile = outfile
        self.params = params["Annotate"]
                
    def callkraken(self):
        kcall=[self.params["Kraken2"]["krakenexe"]]
        kcall.append("--db {}".format(self.params["Kraken2"]["db"]))
        kcall.append("--output {}".format(self.outfile))
        if self.params["Kraken2"]["threads"] != "":
            kcall.append("--threads {}".format(self.params["Kraken2"]["threads"]))
        if self.params["Kraken2"]["additional_args"] != "": 
            kcall.append(self.params["Kraken2"]["additional_args"])
        kcall.append(self.infile)
        return(" ".join(kcall))
        
    
