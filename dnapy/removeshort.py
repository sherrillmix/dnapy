#!/usr/bin/env python
import sys
import argparse
import gzip
import atexit
import Bio.SeqIO.QualityIO
from dnapy import helper




def removeShort(fastqFile,minLength=10):
    nBad=0
    with helper.openNormalOrGz(fastqFile) as fastq:
        print fastq
        for currentRead in Bio.SeqIO.QualityIO.FastqGeneralIterator(fastq):
            if len(currentRead[1])>=minLength:
                yield {"read":currentRead,"nBad":nBad}
            else:
                nBad+=1



def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to remove short reads from a fastq file.")
    parser.add_argument('fastqFile', help='a fastq (potentially gzipped) file containing the alignment',type=helper.check_file)
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument("-l","--minLength", help="minimum length read to output (default:15)",default=15,type=int)
    args=parser.parse_args(argv)
 

    nBad=0
    nGood=0
    for currentRead in removeShort(args.fastqFile,args.minLength):
        helper.writeFastqRead(sys.stdout,currentRead['read'])
        nOld=nBad+nGood
        nBad+=currentRead['nBad']
        nGood+=1
        if args.dots>0:
            nDots=sum([1 for ii in  range(nOld,nGood+nBad+1) if ii%args.dot==0])
            sys.stderr.write('.'*nDots)
        
    if args.dots>0:
        sys.stderr.write("\nGood reads: "+str(nGood)+" Bad reads: "+str(nBad)+"\n")



if __name__ == '__main__':
    main(sys.argv[1:])
