from pyrap.tables import table
import numpy as np
import argparse

def readArguments():
    parser=argparse.ArgumentParser("Script to find empty channels in a given dataset")
    parser.add_argument("-v","--verbose",help="Be verbose, say everything program does. Default is False",required=False,action="store_true")
    parser.add_argument("--filename",type=str,help="Name of the measurement set for which weights want to be calculated",required=True,nargs="+")
    parser.add_argument("--colname",type=str,help="Name of the data column you want to check for emptiness. Default is CORRECTED_DATA.",required=False,default="CORRECTED_DATA")
    parser.add_argument("--flagcol",type=str,help="Name of the flag column you want to apply. Default is FLAG.",required=False,default="FLAG")
    args=parser.parse_args()
    return vars(args)

def ReturnNonEmptyChans(filename,colname="CORRECTED_DATA",flagname="FLAG"):
    # open ms
    ms=table(filename,ack=False)
    # read data, flags
    d=ms.getcol(colname)
    f=ms.getcol(flagname)
    d[f==True]=0.
    # see presence of nonzero vis per chan
    goodchans=[]
    for i in range(d.shape[1]):
        datasum=np.sum(np.abs(d[:,i],dtype=np.float64))
        if datasum!=0.:
            goodchans.append(i)
    return goodchans


### output name dchan_msname

if __name__=="__main__":
    args=readArguments()
    filenames=args["filename"]
    dcol=args["colname"]
    flags=args["flagcol"]
    for filename in filenames:
        goodchans=ReturnNonEmptyChans(filename,colname=dcol,flagname=flags)
        print filename, goodchans
