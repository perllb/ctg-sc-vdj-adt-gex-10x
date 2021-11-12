#!/usr/bin/python3.4
import pandas as pd
import numpy as np
import sys, getopt
import os
import re

def main(argv):
    projdir = ''
    outputdir = ''
    pipeline = ''
    usage='> Usage: ctg-sc-count-mterics-concat.py -i PROJECT-OUTDIR -o SUMMARY-OUTDIR -p PIPELINE'

    try:
        opts, args = getopt.getopt(argv,"hi:o:p:",["projdir=", "outdir=", "pipeline="])
    except getopt.GetoptError:
        print(usage)
        sys.exit(2)
    if len(sys.argv) <= 2:
        print("> Error: No project dir / output dir entered:")
        print(usage)
        sys.exit()
    for opt, arg in opts:
        if opt == '-h':
            print(usage)
            sys.exit()
        elif opt in ("-i", "--projdir"):
            projdir = arg
        elif opt in ("-o", "--outdir"):
            outputdir = arg
        elif opt in ("-p", "--pipeline"):
            pipeline = arg

    print("Outdir: " + outputdir)
    print("Indir: " + projdir)
    projid="CTG-" + os.path.basename(os.path.normpath(projdir))
    print("ProjiD: " + projid)

    # list all metricsfiles
    samples = os.listdir(projdir + '/count/')

    # Iterate through samples, create output csv for each library type
    for sname in samples:

        # Read entire metrics file
        fname = projdir + "/qc/cellranger/" + sname + ".metrics_summary.csv"
        alld = pd.read_csv(fname)

        # Get only Sample rows
        sampd = alld[alld["Library or Sample"]=="Sample"]
        # Fetch all library types of interest
        ltypes=set(list(sampd["Library Type"]))

        # Iterate through all library types, create a file for each
        for lib in ltypes:
            currlib=lib.replace(" ","")
            # Get rows with current lib
            currpd=sampd[sampd["Library Type"]==lib]

            # Prepare outfile for sample and current lib type
            curroutn = outputdir + "/" + sname + "." + currlib + ".metrics_summary.csv"
            currnewpd=pd.DataFrame(columns=["Sample","CTG-Proj-ID"] + list(currpd["Metric Name"]))

            # Add metric values
            currnewpd.loc[0] = [sname,projid] + list(currpd["Metric Value"])

            # Write to csv
            currnewpd.to_csv(curroutn,index=False)

    # For each library type, collect all samples into one table
    for lib in ltypes:
        currlib=lib.replace(" ","")
        out1 = outputdir + "/ctg-cellranger-count-summary_metrics." + currlib + ".csv"
        out2 = outputdir + "/cellranger-count_summary_metrics." + currlib + "_mqc.csv"
        print("Outfile - csv: " + out1)
        print("Outfile - mqc: " + out2)

        # Initiate DF with first sample
        sname = samples[0]
        fname = outputdir + "/" + sname + "." + currlib + ".metrics_summary.csv"
        final = pd.read_csv(fname)

        # Iterate throguh all other samples on current lib and append to table
        for sname in samples[1:]:
            fname = outputdir + "/" + sname + "." + currlib + ".metrics_summary.csv"
            data = pd.read_csv(fname)
            final = pd.concat([final,data],axis=0)
        final=final.set_index("Sample")
        # Remove all commas
        final=final.replace(",","",regex=True)
        final.to_csv(out1,sep=",")

        # Now Samples are appended into on table for current library type - next is to produce the multiqc formatted file
        mqdf = pd.read_csv(out1)

        # Write header
        f = open(out2,'a')
        f.write("# plot_type: 'table'" + "\n")
        f.write("# section_name: '10x Cellranger Metrics: %s'\n" % currlib)
        f.write("# description: '10x %s count metrics'\n" % pipeline)
        f.write("# pconfig:\n")
        f.write("#     namespace: 'CTG'\n")
        f.write("# headers:\n")
        # iterate over column names
        colidx=1
        for col in mqdf:
            f.write("#     col"+str(colidx)+":\n")
            f.write("#         title: '"+col+"'\n")
            f.write("#         description: '"+col+"'\n")
            #check datatype for formating
            form="na"
            if "," in str(mqdf[col][0]):
                f.write("#         format: '{:.0f}' \n")
            if "%" in str(mqdf[col][0]):
                f.write("#         min: 0\n")
                f.write("#         max: 100\n")
                f.write("#         format: '{:.1f}'\n")
                f.write("#         suffix: '%'\n")
            colidx = colidx+1
        f.close()
        # fix colnames
        columns = mqdf.columns
        newcols = ['Sample-ID']
        colidx=2
        for col in range(1,len(columns)):
            newcols.append("col"+str(colidx))
            colidx = colidx+1
        newcols
        mqdf.columns = newcols
        outdf=mqdf.copy(deep=True)
        # Iterate over columns to get format of each
        for rowIndex, row in mqdf.iterrows(): #iterate over rows
            for columnIndex, value in row.items():
                newval=""
                if "%" in str(value):
                    newval=value[:-1]
                    outdf.at[rowIndex,columnIndex] = newval
                if "," in str(value):
                    newval=re.sub(",","",str(value))
                    outdf.at[rowIndex,columnIndex] = newval
        # Write header with format
        outdf.to_csv(out2,mode="a",sep="\t",index=False)


if __name__ == "__main__":
    main(sys.argv[1:])
