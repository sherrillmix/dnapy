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

class splitFastqIter:
    def __init__(self, fastqFile,splits,returnUnassigned=False):
        self.nGood = 0
        self.nBad = 0
        self.fastqFile=fastqFile
        self.splits=splits
        self.fastqHandle=helper.openNormalOrGz(self.fastqFile)
        self.fastq=Bio.SeqIO.QualityIO.FastqGeneralIterator(self.fastqHandle)
        self.returnUnassigned=returnUnassigned
        #make read-file pairs into something efficient for lookup
        self.splitDic=dict(splits)
    
    def __enter__(self):
        return self

    def __exit__(self,exc_type,exc_value,traceback):
        helper.closeFiles(self.fastqHandle)

    def __iter__(self):
        return self

    def next(self):
        return self.__next__()

    def __next__(self):
        for currentRead in self.fastq:
            name=currentRead[0].partition(' ')[0]
            if name in self.splitDic:
                self.nGood+=1
                return (currentRead,self.splitDic[name],True)
            else:
                self.nBad+=1
                if self.returnUnassigned:
                    return (currentRead,None,False)
        raise StopIteration()


def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to take a list of read names and desired files and a fastq file containing various reads and output reads to their assigned file locations. The script takes a standard fastq read file and a csv separated file containing read identifiers and file locations and outputs reads which match the appropriate read names into their assigned files.")
    parser.add_argument('fastqFile', help='a fastq file (potentially gzipped) containing the sequence reads',type=helper.checkFile) #,nargs='+'
    parser.add_argument('-s','--splitFile', help='a (potentially gzipped) file containing comma separated reads names and file names (with no header and no commas in the sample names)',type=helper.checkFile,required=True)
    parser.add_argument("-d","--dots", help="output dot to stderr every X reads. Input a negative number to suppress output (default:-1)", default=-1,type=int)
    parser.add_argument('-o','--outputPath', help='a string giving the desired output directory',type=helper.checkDir,default='.')
    parser.add_argument('-u','--unassigned', help='if set then store unassigned reads to {outputPath}/__UNASSIGNED__.fastq.gz',action='store_true')
    parser.add_argument('-a','--append', help='if set then append to already existing output files',action='store_true')

    args=parser.parse_args(argv)

    splitIds=helper.readSimpleCsv(args.splitFile)
    reads=[xx[0] for xx in splitIds]
    files=[xx[1] for xx in splitIds]
    uniqFiles=set(files)

    if args.append: mode='a'
    else: mode='w'

    #duplicates would be quietly autoreplaced. just let that happen instead of erroring here?
    if(len(set(reads))!=len(reads)):
        raise argparse.ArgumentTypeError("A read appears more than once")
    outputFiles=[os.path.join(args.outputPath,xx) for xx in uniqFiles]
    outHandles=dict(zip(uniqFiles,[helper.openNormalOrGz(xx,mode) for xx in outputFiles]))

    if args.unassigned:
        if '__UNASSIGNED__.fastq.gz' in uniqFiles:
            raise argparse.ArgumentTypeError("Output file named __UNASSIGNED__.fastq.gz clashes with unassigned output. Please rename") #handle better?
        badReadFile=os.path.join(args.outputPath,"__UNASSIGNED__.fastq.gz")
        badReadHandle=helper.openNormalOrGz(badReadFile,mode)

    with splitFastqIter(args.fastqFile,zip(reads,files),args.unassigned) as fastqIter:
        for currentRead,splitId,assigned in fastqIter:
            if args.unassigned and not assigned:
                helper.writeFastqRead(badReadHandle,currentRead)
            else:
                helper.writeFastqRead(outHandles[splitId],currentRead)
                if args.dots>0:
                    if fastqIter.nGood % args.dots==0:
                        sys.stderr.write('.')
                        sys.stderr.flush()

        if args.dots>0:
            sys.stderr.write("\nReads assigned to files: "+str(fastqIter.nGood)+" Unassigned reads: "+str(fastqIter.nBad)+"\n")

    helper.closeFiles(outHandles)
    #for key in outHandles:
        #outHandles[key].close()


if __name__ == '__main__':
    main()
