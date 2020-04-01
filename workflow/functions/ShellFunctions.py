#general shell functions
from workflow.functions import SequencingData as SD
import os

#function to symlink input file to output
def symlink(inf,outf):
    f=SD.SequencingData(inf)
    command="ln -rs {} {}".format(inf,outf)
    if f.pairedname != "":
        outp=os.path.dirname(outf)+"/"+f.pairedname
        command+=" && ln -rs {} {}".format(f.pairedname,outp)
    return(command)
