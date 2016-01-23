#!/usr/bin/env python
import sys
import argparse
import gzip
import atexit
import Bio.SeqIO.QualityIO
from dnapy import helper






class shortFilterFastqIter:
    def __init__(self, fastqFile, minLength=10):
        self.nGood = 0
        self.nBad = 0
        self.minLength=minLength
        self.fastqFile=fastqFile
        self.fastqHandle=helper.openNormalOrGz(self.fastqFile)
        self.fastq=Bio.SeqIO.QualityIO.FastqGeneralIterator(self.fastqHandle)
    
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
                self.nGood+=1
                return currentRead
            else:
                self.nBad+=1
        raise StopIteration()


#def removeShort(fastqFile,minLength=10):
    #nBad=0
    #fastqIter=shortFilterFastqIter(fastqFile)
    #for currentRead in fastqIter:
        #
    #with helper.openNormalOrGz(fastqFile) as fastq:
        #for currentRead in Bio.SeqIO.QualityIO.FastqGeneralIterator(fastq):
            #if len(currentRead[1])>=minLength:
                #yield {"read":currentRead,"nBad":nBad}
            #else:
                #nBad+=1



def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to remove short reads from a fastq file.")
    parser.add_argument('fastqFile', help='a fastq (potentially gzipped) file containing the alignment',type=helper.check_file)
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument("-l","--minLength", help="minimum length read to output (default:15)",default=15,type=int)
    args=parser.parse_args(argv)
 

    with shortFilterFastqIter(args.fastqFile,args.minLength) as fastqIter:
        for currentRead in fastqIter:
            helper.writeFastqRead(sys.stdout,currentRead)
            if args.dots>0:
                if fastqIter.nGood % args.dots==0:
                    sys.stderr.write('.')

        if args.dots>0:
            sys.stderr.write("\nGood reads: "+str(fastqIter.nGood)+" Bad reads: "+str(fastqIter.nBad)+"\n")



if __name__ == '__main__':
    main(sys.argv[1:])
