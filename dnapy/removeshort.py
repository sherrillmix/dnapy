#!/usr/bin/env python
import sys
import argparse
import gzip
import atexit
import Bio.SeqIO.QualityIO
from dnapy import helper



def closeFiles(openFiles):
    for openFile in openFiles:
        openFile.close()

def writeRead(readFile,read):
    readFile.write("@%s\n%s\n+%s\n%s\n" %(read[0],read[1],read[0],read[2]))

def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to remove short reads from a fastq file.")
    parser.add_argument('fastqFile', help='a fastq (potentially gzipped) file containing the alignment',type=helper.check_file)
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument("-l","--minLength", help="minimum length read to output (default:15)",default=15,type=int)
    args=parser.parse_args(argv)
 

    if len(args)!=1:
        sys.stderr.write('Please provide a single fastq or fastq.gz file\n')
        return 4
    try:
        if args[0][-2:]=='gz' or args[0][-4]=='gzip':
            fastq=gzip.open(args[0], "r")
        else:
            fastq=open(args[0],"r")
    except IOError:
        sys.stderr.write("Problem opening file:"+args[0]+"\n")
        return 3
    openFiles=[fastq]
    atexit.register(closeFiles,openFiles)
    

    for currentRead in Bio.SeqIO.QualityIO.FastqGeneralIterator(fastq):
        if len(currentRead[1])>=args.minLength:
            nGood+=1
            writeRead(sys.stdout,currentRead)
        else:
            nBad+=1
        if args.dots>0 and (nGood+nBad) % args.dots == 0:
            sys.stderr.write('.')
    if verbose:
        sys.stderr.write("\nGood reads: "+str(nGood)+" Bad reads: "+str(nBad)+"\n")


if __name__ == '__main__':
    main(sys.argv[1:])
