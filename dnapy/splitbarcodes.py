
import sys
import argparse
import gzip
import Bio.SeqIO.QualityIO
from dnapy import helper
import os

#make sure zip is iteratable
try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass

class barcodeFastqIter:
    def __init__(self, fastqFiles,indexFiles,barcodes):
        self.nGood = 0
        self.nBad = 0
        self.fastqFiles=fastqFiles
        self.indexFiles=indexFiles
        self.barcodes=barcodes
        self.fastqHandles=[helper.openNormalOrGz(x) for x in self.fastqFiles]
        self.indexHandles=[helper.openNormalOrGz(x) for x in self.indexFiles]
        self.fastqs=[Bio.SeqIO.QualityIO.FastqGeneralIterator(x) for x in self.fastqHandles]
        self.indexs=[Bio.SeqIO.QualityIO.FastqGeneralIterator(x) for x in self.indexHandles]
    
    def __enter__(self):
        return self

    def __exit__(self,exc_type,exc_value,traceback):
        helper.closeFiles(self.fastqHandles)

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        for currentReads,currentBars in zip(zip(*self.fastqs),zip(*self.indexs)):
            bars=tuple([x[1] for x in currentBars])
            if bars in self.barcodes:
                self.nGood+=1
                return (currentReads,bars)
            else:
                self.nBad+=1
        raise StopIteration()


def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to take a list of barcodes and one or more fastq reads and one or two index reads and output reads matching the barcodes into a seperate file for each barcode. The script takes read files and index files where the reads and indexs are in the same order and outputs reads which match the appropriate barcodes into separate files.")
    parser.add_argument('fastqFiles', help='a fastq file(s) (potentially gzipped) containing the sequence reads',type=helper.checkFile,nargs='+')
    parser.add_argument('-i','--indexFiles', help='a fastq file(s) (potentially gzipped) containing the index reads',type=helper.checkFile,nargs='+')
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument('-b','--barcodeFile', help='a (potentially gzipped) file containing comma separated sample names, first barcode and second barcode (with no header and no commas in the sample names)',type=helper.checkFile,required=True)
    parser.add_argument('-o','--outputPath', help='a string giving the desired output directory',type=helper.checkDir,default='.')

    args=parser.parse_args(argv)

    nFiles=len(args.fastqFiles)
    barcodes=helper.readSimpleCsv(args.barcodeFile)
    if(len(args.indexFiles)!=len(barcodes[0])-1):
        raise argparse.ArgumentTypeError("Number of index files and index columns in the barcodeFile do not agree")

    samples=[xx[0] for xx in barcodes]
    if(len(set(samples))!=len(samples)):
        raise argparse.ArgumentTypeError("Two or more samples share the same name")
    outputFiles=[[os.path.join(args.outputPath,xx)+"_"+str(ii+1)+".fastq.gz" for ii in range(nFiles)] for xx in samples]
    bars=[tuple(xx[1:]) for xx in barcodes]
    barSet=set(bars)
    if len(barSet)!=len(bars):
        raise argparse.ArgumentTypeError("Two or more samples share the same set of barcodes")
    outHandles=dict(zip(bars,[[helper.openNormalOrGz(yy,'w') for yy in xx] for xx in outputFiles]))


    with barcodeFastqIter(args.fastqFiles,args.indexFiles,bars) as fastqIter:
        for currentReads,bar in fastqIter:
            for read,outFile in zip(currentReads,outHandles[bar]):
                helper.writeFastqRead(outFile,read)
            if args.dots>0:
                if fastqIter.nGood % args.dots==0:
                    sys.stderr.write('.')
                    sys.stderr.flush()

        if args.dots>0:
            sys.stderr.write("\nReads assigned to barcode: "+str(fastqIter.nGood)+" Unassigned reads: "+str(fastqIter.nBad)+"\n")

    for key in outHandles:
        helper.closeFiles(outHandles[key])


if __name__ == '__main__':
    main(sys.argv[1:])
