#!/usr/bin/env python
import sys
import argparse
import gzip
import atexit
import Bio.SeqIO.QualityIO
from dnapy import helper
import re



class shortFilterFastqIter:
    def __init__(self, fastqFile, minLength=10, removeN=False, suppressBad=False):
        self.nGood = 0
        self.nBad = 0
        self.minLength=minLength
        self.removeN=removeN
        self.suppressBad=suppressBad
        self.fastqFile=fastqFile
        self.fastqHandle=helper.openNormalOrGz(self.fastqFile)
        if suppressBad:
            self.fastq=helper.readSimpleFastq(self.fastqHandle)
        else:
            self.fastq=Bio.SeqIO.QualityIO.FastqGeneralIterator(self.fastqHandle)
        self.nSearch=re.compile('[^ACTG]').search
        self.badList=[]
    
    def __enter__(self):
        return self

    def __exit__(self,exc_type,exc_value,traceback):
        self.fastqHandle.close()

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        for currentRead in self.fastq:
            if len(currentRead[1])>=self.minLength:
                if not self.removeN or not bool(self.nSearch(currentRead[1])):
                    if not self.suppressBad or len(currentRead[1])==len(currentRead[2]):
                        self.nGood+=1
                        return currentRead
                    else:
                        if len(self.badList)<10000:
                            self.badList.append(currentRead)
            self.nBad+=1
        raise StopIteration()


def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to remove short reads from a fastq file.")
    parser.add_argument('fastqFile', help='a (potentially gzipped) fastq file containing the sequence data',type=helper.checkFile)
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument("-l","--minLength", help="minimum length read to output (default:15)",default=15,type=int)
    parser.add_argument("-n","--removeN", help="remove reads which contain anything other than A, C, T or G",action='store_true')
    parser.add_argument("-p","--removePoor", help="remove reads with different length sequence and qualities.  Note this requires assuming that all reads are 4 lines each",action='store_true')
    parser.add_argument("-b","--badOut", help="a file path in which to save the first 10000 malformed reads with different length sequence and qualities. Note this requires assuming that all reads are 4 lines each",type=argparse.FileType('w'),default=None)
    args=parser.parse_args(argv)
 

    with shortFilterFastqIter(args.fastqFile,args.minLength,args.removeN,args.removePoor) as fastqIter:
        for currentRead in fastqIter:
            helper.writeFastqRead(sys.stdout,currentRead)
            if args.dots>0:
                if fastqIter.nGood % args.dots==0:
                    sys.stderr.write('.')
                    sys.stderr.flush()

        if args.badOut is not None:
            for currentRead in fastqIter.badList:
                helper.writeFastqRead(args.badOut,currentRead)
            args.badOut.close()

        if args.dots>0:
            sys.stderr.write("\nGood reads: "+str(fastqIter.nGood)+" Bad reads: "+str(fastqIter.nBad)+"\n")



if __name__ == '__main__':
    main()
