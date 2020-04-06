#function to make a bash call to count the number of reads in each filtered file
from workflow.functions import SequencingData as SD

def CountFiltered(orig,rrna,host,outfile,params):
    #generate awk to count reads
    def counter(seqfile,outfile):
        sdat=SD.SequencingData(seqfile)
        div=4
        if sdat.fileformat == "fasta":
            div=2
        return("""zcat {0} | wc -l | awk 'BEGIN {1}ORS=\"\"{2}; END {1}x=$1/{3}; print \"\\t\"x{2}' >> {4}""".format(seqfile,"{{","}}",div,outfile))
    #get counts at start, after rrna and after host filtering
    original = SD.SequencingData(orig)
    call=["""printf "File\\tInput\\tPost_rRNA_Filtering\\tPost_Genome_Filtering\\n{}" > {}""".format(original.cleanname,outfile)]
    ocount = counter(orig,outfile)
    rcount = """printf "\\tNA" >> {}""".format(outfile)
    gcount = """printf "\\tNA" >> {}""".format(outfile)
    if params["Filter"]["General"]["rrna_filter"] == "true": 
        rcount = counter(rrna,outfile)
    if params["Filter"]["General"]["host_filter"] == "true":
        gcount = counter(host,outfile)
    call.append(ocount)
    call.append(rcount)
    call.append(gcount)
    call.append("""printf "\\n" >> {}""".format(outfile))
    return(""" && """.join(call))
