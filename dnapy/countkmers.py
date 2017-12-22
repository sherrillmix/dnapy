#!/usr/bin/env python
import collections
from dnapy import helper
import argparse
import Bio.SeqIO.QualityIO
import sys
import multiprocessing

#make sure zip is iteratable
try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass

def splitIntoKmer(read,k=10,step=None):
    if step is None: step=k
    starts=range(0,len(read)-k+1,step)
    return([read[start:start+k] for start in starts])

def generateAllKmers(k=10,bases=['A','C','G','T']):
    if k < 1: return([])
    if k == 1: return bases
    suffixs=generateAllKmers(k=k-1,bases=bases)
    return([base+suffix for base in bases for suffix in suffixs])

def countKmersInReads(reads,k=10):
    kmers=collections.defaultdict(defaultVal)
    for read in reads:
        readKmers=splitIntoKmer(read[1],k)
        for kmer in readKmers:
            kmers[kmer]+=1
    return(kmers)

def countKmersInFile(args):
    fastqFile=args[0]
    k=args[1]
    sys.stderr.write('Working on file %s\n' % fastqFile)
    with helper.openNormalOrGz(fastqFile) as fastqHandle:
        fastq=Bio.SeqIO.QualityIO.FastqGeneralIterator(fastqHandle)
        return(countKmersInReads(fastq,k))

#can't use lambda or break pickel and multiprocess
def defaultVal():
    return 0

def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to take a fastq file(s) and count the total k-mers across all reads in each file. Note that partial kmers are discarded e.g. the last 3 reads of a 23 base read will be ignored. Return a comma separated file with a header row then a row for each file and a file column then a column for each kmer")
    parser.add_argument('fastqFiles', help='a fastq file(s) (potentially gzipped) containing the sequence reads',type=helper.checkFile,nargs='+')
    parser.add_argument('-k','--kmerLength', help='the lengh of kmer to be used. Be careful with values larger than 20.',default=10,type=helper.checkPositiveInt)
    parser.add_argument('-t','--nThreads', help='the number of threadss to use for processing. Should be less than or equal the number of threads on computer.',default=1,type=helper.checkPositiveInt)

    args=parser.parse_args(argv)
    nFiles=len(args.fastqFiles)

    #kmerCounts=[[] for _ in range(nFiles)]
    p=multiprocessing.Pool(args.nThreads)
    kmerCounts=p.map(countKmersInFile,zip(args.fastqFiles,[args.kmerLength]*nFiles))

    presentKmers=sorted(set(sum([list(xx.keys()) for xx in kmerCounts],[])))

    #allKmers=generateAllKmers(args.kmerLength)
    #print("file,%s" % ",".join(allKmers))
    print("kmer,%s" % ",".join(args.fastqFiles))
    for kmer in presentKmers:
        print("%s,%s" % (kmer,",".join([str(kmerCount[kmer]) for kmerCount in kmerCounts])))

if __name__ == '__main__':
    main(sys.argv[1:])
