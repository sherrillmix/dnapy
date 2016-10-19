
import sys
import argparse
import gzip
import Bio.SeqIO.QualityIO
from dnapy import helper


def readBarcodes(barFile):
    lines=line.strip() for line in helper.openNormalOrGz(args.filterFile)
    #TODO SPLIT BY COMMA
    # return dictionary

class barodeFastqIter:
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
            names=[x[0].partition(' ')[0] for x in currentReads]
            bars=[x[1] for x in currentBars]
            #TODO ASSIGN BARCODE
            if not any([x in self.patterns for x in names]):
                self.nGood+=1
                return currentReads
            else:
                self.nBad+=1
        raise StopIteration()

def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to take a list of barcodes and one or more fastq reads and one or two index reads and output reads matching the barcodes into a seperate file for each barcode. The script looks for reads which have a name line where the string before a space exactly matches a pattern. If multiple files are passed in, then they are processed in sync and if any name matches that read is discarded from all files.")
    parser.add_argument('fastqFiles', help='a fastq or fastqs (potentially gzipped) file containing the alignment',type=helper.check_file,nargs='+')
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument('-b','--barcodeFile', help='a file (potentially gzipped) file containing comma separated sample names, first barcode and second barcode',type=helper.check_file,required=True)
    #TODO CHECK IF REAL DIRECTORY/CREATE?
    parser.add_argument('-o','--outputPath', help='a string giving the desired output directory. (default='.')',type=str)

    args=parser.parse_args(argv)
    if(args.outputFiles is None):
        outputFiles=['out'+str(ii)+'.fastq.gz' for ii in range(1,len(args.fastqFiles)+1)]
    else:
        outputFiles=args.outputFiles.split(',')
    if(len(outputFiles)!=len(args.fastqFiles)):
        argparse.ArgumentTypeError("Input and output file numbers do not match")
    outHandles=[helper.openNormalOrGz(x,'w') for x in outputFiles]

    barcodes1=readBarcodes(args.barcodeFile)

    with filterFastqIter(args.fastqFiles,patterns) as fastqIter:
        for currentReads in fastqIter:
            for read,outFile in zip(currentReads,outHandles):
                helper.writeFastqRead(outFile,read)
            if args.dots>0:
                if fastqIter.nGood % args.dots==0:
                    sys.stderr.write('.')

        if args.dots>0:
            sys.stderr.write("\nGood reads: "+str(fastqIter.nGood)+" Bad reads: "+str(fastqIter.nBad)+"\n")

    helper.closeFiles(outHandles)


if __name__ == '__main__':
    main(sys.argv[1:])
