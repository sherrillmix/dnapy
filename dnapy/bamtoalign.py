import argparse
import pysam
import sys
from dnapy import helper



def getAlignsInFile(inputFile,region=None):
    with pysam.AlignmentFile(inputFile, "rb" ) as samfile:
        for read in samfile.fetch(region=region):
            rPos=read.reference_start
            qPos=0
            for operation,length in read.cigartuples:
                print operation
                print length
                if operation in [0,1,4,7,8]:
                    qPos+=length
                    print qPos
            print read.query_sequence
            print read.get_reference_positions()
            #filter here
            yield {"name": read.query_name, "start": read.pos+1, "end": read.aend,'strand': '-' if read.flag&16==16 else '+'}



def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to convert a bam file into an aligned fasta file. The command generates fasta formatted output (two lines for each sequence: a name line prepended by > and a line containing the aligned sequence) to standard out.")
    parser.add_argument('bamFile', help='a bam file containing the alignment',type=helper.checkFile)
    parser.add_argument("-v","--verbose", help="increase output verbosity to stderr", action="store_true")
    parser.add_argument("-r","--region", help="the region to pull reads from",default=None)
    parser.add_argument("-e","--endSpan", help="ignore spans of matches at the start or end of a read less than this cutoff",default=0,type=int)
    args=parser.parse_args(argv)
 
        
    if args.verbose:
        sys.stderr.write("Arguments: ")
        for key, value in vars(args).items():
            sys.stderr.write("   "+key+": "+str(value))

    for read in getAlignsInFile(args.bamFile,args.region):
        print read
        break
            

    
if __name__ == '__main__':
    main()
