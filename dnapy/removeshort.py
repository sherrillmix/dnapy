#!/usr/bin/env python
import sys
import argparse
import gzip
import atexit
import Bio.SeqIO.QualityIO
from dnapy import helper




def removeShort(fastqFile,minLength=10,dots=0):
    nGood=0
    nBad=0
    with helper.openNormalOrGz(fastqFile) as fastq:
        for currentRead in Bio.SeqIO.QualityIO.FastqGeneralIterator(fastq):
            if len(currentRead[1])>=minLength:
                nGood+=1
                helper.writeFastqRead(sys.stdout,currentRead)
            else:
                nBad+=1
            if dots>0 and (nGood+nBad) % dots == 0:
                sys.stderr.write('.')
        if dots>0:
            sys.stderr.write("\nGood reads: "+str(nGood)+" Bad reads: "+str(nBad)+"\n")



def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to remove short reads from a fastq file.")
    parser.add_argument('fastqFile', help='a fastq (potentially gzipped) file containing the alignment',type=helper.check_file)
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument("-l","--minLength", help="minimum length read to output (default:15)",default=15,type=int)
    args=parser.parse_args(argv)
 

    removeShort(fastq,args.minLength,args.dots)



if __name__ == '__main__':
    main(sys.argv[1:])
