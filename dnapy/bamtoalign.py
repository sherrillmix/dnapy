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
    parser = argparse.ArgumentParser(description="A program to convert a bam file into an aligned fasta file. The command generates fasta formatted output (two lines for each sequence: a name line prepended by > and a line containing the aligned sequence) to standard out.")
    parser.add_argument('bamFile', help='a bam file containing the alignment',type=helper.checkFile)
    parser.add_argument("-v","--verbose", help="increase output verbosity to stderr", action="store_true")
    parser.add_argument("-r","--region", help="the region to pull reads from",default=None)
    parser.add_argument("-e","--endSpan", help="Ignore spans of matches at the start or end of a read less than this cutoff",default=0,type=int)
    args=parser.parse_args(argv)
 
        
    if args.verbose:
        sys.stderr.write("Arguments: ")
        for key, value in vars(args).items():
            sys.stderr.write("   "+key+": "+str(value))

    if args.file is None:
        regions=[args.region]
    else:
        with open(args.file) as f:
            regions=[x.strip('\n') for x in f.readlines()]


    if(not args.noHeader):
        if args.regionColumn: header="region,ref,start,end,strand"
        else: header="ref,start,end,strand"
        print(header)
    for region in regions:
        for read in getStartsInFile(args.bamFile,region,args.maxGaps):
            if args.regionColumn:
                print("%s,%s,%d,%d,%s" % (region,read["ref"],read['start'],read['end'],read['strand']))
            else:
                print("%s,%d,%d,%s" % (read["ref"],read['start'],read['end'],read['strand']))
            

    
if __name__ == '__main__':
    main()
