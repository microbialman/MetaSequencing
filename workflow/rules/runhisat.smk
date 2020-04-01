#Rule to run HISAT2

#import class to generate the call to hisat
from workflow.functions import  as F

rule runhisat:
    input:
        
