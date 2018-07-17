#!/usr/bin/env python
import sys
import argparse
import gzip
import Bio.SeqIO.QualityIO
from dnapy import helper

#make sure zip is iteratable
try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass

class filterFastqIter:
    def __init__(self, fastqFiles,patterns,keep=False):
        self.nGood = 0
        self.nBad = 0
        self.fastqFiles=fastqFiles
        self.patterns=patterns
        self.fastqHandles=[helper.openNormalOrGz(x) for x in self.fastqFiles]
        self.fastqs=[Bio.SeqIO.QualityIO.FastqGeneralIterator(x) for x in self.fastqHandles]
        self.keep=keep
    
    def __enter__(self):
        return self

    def __exit__(self,exc_type,exc_value,traceback):
        helper.closeFiles(self.fastqHandles)

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        for currentReads in zip(*self.fastqs):
            names=[x[0].partition(' ')[0] for x in currentReads]
            if self.keep ^ (not any([x in self.patterns for x in names])):
                self.nGood+=1
                return currentReads
            else:
                self.nBad+=1
        raise StopIteration()


def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to filter reads by name from a single/set of fastq file(s). The script looks for reads which have a name line where the string before a space exactly matches a pattern. If multiple files are passed in, then they are processed in sync and if any name matches that read is discarded (or kept) from all files.")
    parser.add_argument('fastqFiles', help='a (potentially gzipped) fastq file(s) containing the reads with the order of reads the same in all files',type=helper.checkFile,nargs='+')
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument('-o','--outputFiles', help='an output file(s) (one for each input fastq file). default(out1.fastq.gz ... outn.fastq.gz where n is the number of fastqFiles)',type=str,nargs='*')
    parser.add_argument('-f','--filterFile', help='a (potentially gzipped) file containing the names of reads to be filtered one per line',type=helper.checkFile,required=True)
    parser.add_argument('-k','--keep', help='keep reads matching the filter file and filter all nonmatching reads',action='store_true')


    args=parser.parse_args(argv)
    if(args.outputFiles is None):
        outputFiles=['out'+str(ii)+'.fastq.gz' for ii in range(1,len(args.fastqFiles)+1)]
    else:
        outputFiles=args.outputFiles
    if(len(outputFiles)!=len(args.fastqFiles)):
        raise argparse.ArgumentTypeError("Input and output file numbers do not match")
    outHandles=[helper.openNormalOrGz(x,'w') for x in outputFiles]

    patterns=set(line.strip() for line in helper.openNormalOrGz(args.filterFile))

    with filterFastqIter(args.fastqFiles,patterns,args.keep) as fastqIter:
        for currentReads in fastqIter:
            for read,outFile in zip(currentReads,outHandles):
                helper.writeFastqRead(outFile,read)
            if args.dots>0:
                if fastqIter.nGood % args.dots==0:
                    sys.stderr.write('.')
                    sys.stderr.flush()

        if args.dots>0:
            sys.stderr.write("\nGood reads: "+str(fastqIter.nGood)+" Bad reads: "+str(fastqIter.nBad)+"\n")

    helper.closeFiles(outHandles)


if __name__ == '__main__':
    main()
