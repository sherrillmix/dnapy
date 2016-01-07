import argparse
import pysam
import collections
import os
import sys
from dnapy import helper


MAX_DEPTH=1e9

def countBasesInFile(inputFile,region=None):
    with pysam.AlignmentFile(inputFile, "rb" ) as samfile:
        #print(samfile.pileup(region=region))
        for pileupcolumn in samfile.pileup(region=region,max_depth=MAX_DEPTH):
            #print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
            counts=collections.defaultdict(lambda:0)
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    thisBase=pileupread.alignment.query_sequence[pileupread.query_position]
                    thisStrand='-' if pileupread.is_reverse else '+'
                    #counts[thisBase+thisStrand]+=1
                    counts[thisBase]+=1
            yield {"ref": pileupcolumn.reference_name, "pos": pileupcolumn.pos, "n": pileupcolumn.n, "A": counts['A'], "C": counts['C'], "G": counts['G'], "T":counts['T']}



def main(argv=None):
    if argv is None:
        argv=sys.argv[1:]

    parser = argparse.ArgumentParser(description="A program to count the number of bases at each position in a region. The command generates standard output with columns referenceName, position, numberOfReads, and numbers of A, C, G, T.")
    parser.add_argument("-v","--verbose", help="increase output verbosity to stderr", action="store_true")
    parser.add_argument("-b","--bamFile", help="a bam file containing the alignment",type=helper.check_file,required=True)
    parser.add_argument("-r","--region", help="the region to count in",default=None)
    args=parser.parse_args(argv)
        
    if args.verbose:
        sys.stderr.write("Arguments: ")
        for key, value in vars(args).items():
            sys.stderr.write("   "+key+": "+str(value))

        for column in countBasesInFile(args.bamFile,args.region):
            print("%s,%d,%d,%d,%d,%d,%d" % (column["ref"],column['pos'],column['n'],column["A"],column['C'],column['G'],column["T"]))
            

    
if __name__ == '__main__':
    main(sys.argv[1:])
