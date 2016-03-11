import argparse
import pysam
import sys
from dnapy import helper



def getStartsInFile(inputFile,region=None,maxGaps=0):
    with pysam.AlignmentFile(inputFile, "rb" ) as samfile:
        #print(samfile.pileup(region=region))
        for read in samfile.fetch(region=region):
            #filter here
            gaps=sum([length for op, length in read.cigar if op in [1,2,3]]) #insertion,deletion,refskip
            if gaps<=maxGaps:
                yield {"ref": read.reference_name, "start": read.pos+1, "end": read.aend,'strand': '-' if read.flag&16==16 else '+'}



def main(argv=None):
    if argv is None:
        argv=sys.argv[1:]

    parser = argparse.ArgumentParser(description="A program to pull start and end positions in a region. The command generates standard output with columns referenceName, start (1-based), end (1-based), strand ")
    parser.add_argument('bamFile', help='a bam file containing the alignment',type=helper.check_file)
    parser.add_argument("-v","--verbose", help="increase output verbosity to stderr", action="store_true")
    parser.add_argument("-g","--maxGaps", help="maximum allowed insertions or deletions in a read. Otherwise discard", type=int,default=0)
    parser.add_argument("-r","--region", help="the region to count in",default=None)
    args=parser.parse_args(argv)
        
    if args.verbose:
        sys.stderr.write("Arguments: ")
        for key, value in vars(args).items():
            sys.stderr.write("   "+key+": "+str(value))

    if args.region not None and args.file not None:
        raise argparse.ArgumentError("Please specify either a region or region file but not both")

    if args.file is None:
        regions=[args.region]


    header="ref,start,end,strand"
    print(header)
    for region in regions:
        for read in getStartsInFile(args.bamFile,region,args.maxGaps):
            if args.firstCol:
                print("%s,%s,%d,%d,%s" % (region,read["ref"],read['start'],read['end'],read['strand']))
            else:
                print("%s,%d,%d,%s" % (read["ref"],read['start'],read['end'],read['strand']))
            

    
if __name__ == '__main__':
    main(sys.argv[1:])
