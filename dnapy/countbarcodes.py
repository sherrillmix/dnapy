import sys
import argparse
import Bio.SeqIO.QualityIO
from dnapy import helper


class barcodeFastqIter:
    def __init__(self, fastqFiles,start,end,barcodes):
        self.nGood = 0
        self.nBad = 0
        self.fastqFiles=fastqFiles
        self.start=start
        self.end=end
        self.onlyInSet=len(barcodes)>0
        self.barcodes=set(barcodes)
        self.fastqHandles=[helper.openNormalOrGz(x) for x in self.fastqFiles]
        self.fastqs=[Bio.SeqIO.QualityIO.FastqGeneralIterator(x) for x in self.fastqHandles]
    
    def __enter__(self):
        return self

    def __exit__(self,exc_type,exc_value,traceback):
        helper.closeFiles(self.fastqHandles)

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        for fastq in self.fastqs:
            for currentRead in fastq:
                barcode=currentRead[1][(self.start-1):(self.end)]#cut out barcode[x[0].partition(' ')[0] for x in currentReads]
                if not self.onlyInSet or barcode in self.barcodes:
                    self.nGood+=1
                    return barcode
                else:
                    self.nBad+=1
        raise StopIteration()

def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to take a fastq file and count barcodes (potentially only including those found in a whitelist). The script takes a standard fastq read file and an optional barcode whitelist file containing one barcode per line and outputs a header-less csv with columns barcode, count for each barcode to standard out")
    parser.add_argument('fastqFiles', help='a fastq file(s) (potentially gzipped) containing the sequence reads',type=helper.checkFile,nargs='+')
    parser.add_argument("-s","--start", help="Start position of barcode in reads (1-based)", default=1,type=int)
    parser.add_argument("-e","--end", help="End position of barcode in reads (1-based)", default=16,type=int)
    parser.add_argument('-w','--whitelist', help='a (potentially gzipped) file containing whitelisted barcodes one to a line with no header',type=helper.checkFile)
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument("-a","--all", help="if set output 0 for all barcodes in whitelist not appearing", action='store_true')


    args=parser.parse_args(argv)


    if args.whitelist:
        barcodes=[xx[0] for xx in helper.readSimpleCsv(args.whitelist)]
    else:
        barcodes=[]
    #error checking?
    counts={}
    if args.all:
        counts={xx:0 for xx in barcodes}
    with barcodeFastqIter(args.fastqFiles,args.start,args.end,barcodes) as fastqIter:
        for barcode in fastqIter:
            if barcode in counts:
                counts[barcode]+=1
            else:
                counts[barcode]=1
            if args.dots>0:
                if fastqIter.nGood % args.dots==0:
                    sys.stderr.write('.')
                    sys.stderr.flush()
    if args.dots>0:
        sys.stderr.write("\nMatching barcodes found: "+str(fastqIter.nGood)+" Unassigned barcodes: "+str(fastqIter.nBad)+"\n")
    for barcode,count in counts.items():
        sys.stdout.write('%s,%d\n' % (barcode,count))


if __name__ == '__main__':
    main()
