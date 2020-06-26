'''
classes to build calls to Eggnogmapper
'''

class Eggnogmapper:    
    def __init__(self,infile,outfile,params,call):
        self.infile =infile
        self.outfile = outfile
        self.params = params["Annotate"]
        self.statementlist = []
        self.executable="emapper.py"
        self.callemapper()
        if call == "ortho":
            self.orthocall()
        if call == "annot":
            self.annotcall()
            
    def callemapper(self):
       if self.params["Eggnogmapper"]["preload"] != "":
           self.statementlist.append(self.params["Eggnogmapper"]["preload"])
       if self.params["Eggnogmapper"]["eggexe"] != "":
           self.executable=self.params["Eggnogmapper"]["eggexe"]
           
    def orthocall(self):
        oc=[]
        oc.append(self.executable)
        oc.append("-i {} -o {}".format(self.infile,self.outfile.replace(".emapper.seed_orthologs","")))
        oc.append("--data_dir {}".format(self.params["Eggnogmapper"]["eggdata"]))
        oc.append("-m diamond --no_annot --no_file_comments")
        oc.append("--cpu {}".format(self.params["Eggnogmapper"]["threads_orth"]))
        if self.params["Eggnogmapper"]["additional_args_orth"] != "":
            oc.append(self.params["Eggnogmapper"]["additional_args_orth"])        
        self.statementlist.append(" ".join(oc))

    def annotcall(self):
        ac=[]
        ac.append(self.executable)
        ac.append("--annotate_hits_table {} -o {} --no_file_comments".format(self.infile,self.outfile.replace(".emapper.annotations","")))
        ac.append("--cpu {} --data_dir {}".format(self.params["Eggnogmapper"]["threads_annot"],self.params["Eggnogmapper"]["eggdata"]))
        self.statementlist.append(" ".join(ac))

    def build(self):
        return(" && ".join(self.statementlist))
        
