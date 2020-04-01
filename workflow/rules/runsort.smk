#function to generate the call to sortmerna
from workflow.functions import SortMeRNA as S

rule runsortmerna:
    input:
        "{file}"
    output:
        "rrna_filter_out.dir/non_rrna/{file}"
    threads:
        int(config["SortMeRNA"]["threads"])
    resources:
        mem=int(config["SortMeRNA"]["memory"])*1000
    run:
        for i in range(len(input)):
            if config["General"]["rrna_filter"] == "true":
                command=S.SortMeRNA(input[i],output[i],config).build()
            else:
                command="ln -s {} {}".format(input[i],output[i])
            shell(command)
