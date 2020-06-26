#general shell functions
from workflow.functions import SequencingData as SD
import os
import re

#function to symlink input file to output
def symlink(inf,outf):
    f=SD.SequencingData(inf)
    command="ln -rs {} {}".format(inf,outf)
    if f.pairedname != "":
        outp=os.path.dirname(outf)+"/"+f.pairedname
        command+=" && ln -rs {} {}".format(f.pairedname,outp)
    return(command)

#function to generate call to merge contig summary files
def mergecontigsumm(inputs,output,params):
    statementlist = []
    in0 = inputs[0] 
    statementlist.append("> {}".format(output))
    statementlist.append("head -1 {} >>{}".format(in0,output))
    statementlist.append("sed  -i '1s/^/{}\\t{}\\t /' {}".format("file","assembler",output))
    assem=params["Assemble"]["General"]["assembler"]
    #extract filenames to add to summary text file
    for infile in inputs:
        filen=re.search("contig_summaries/(.*).contigs.summary",infile).group(1)
        #just append the last line and add filename and assembler name
        statementlist.append("tail -1 {} >> {}".format(infile,output))
        statementlist.append("sed -i '$s/^/{}\\t{}\\t /' {}".format(filen,assem,output))
    statement = " && ".join(statementlist)
    return(statement)
