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
            counts={'+':collections.defaultdict(lambda:0),'-':collections.defaultdict(lambda:0)}
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    thisBase=pileupread.alignment.query_sequence[pileupread.query_position]
                    thisStrand='-' if pileupread.alignment.is_reverse else '+'
                    counts[thisStrand][thisBase]+=1
                    #counts[thisBase]+=1
            yield {"ref": pileupcolumn.reference_name, "pos": pileupcolumn.pos, "n": pileupcolumn.n, "+": counts['+'], "-": counts['-']}



def main(argv=None):
    bases=['A','C','G','T']
    strands=['+','-']
    if argv is None:
        argv=sys.argv[1:]

    parser = argparse.ArgumentParser(description="A program to count the number of bases at each position in a region. The command generates standard output with columns referenceName, position, numberOfReads, and numbers of A, C, G, T (or A+, A-, C+, ... if --strand).")
    parser.add_argument('bamFile', help='a bam file containing the alignment',type=helper.check_file)
    parser.add_argument("-v","--verbose", help="increase output verbosity to stderr", action="store_true")
    parser.add_argument("-r","--region", help="the region to count in",default=None)
    parser.add_argument("-s","--strand", help="break base counts into positive and negative strand alignments", action="store_true")
    args=parser.parse_args(argv)
        
    if args.verbose:
        sys.stderr.write("Arguments: ")
        for key, value in vars(args).items():
            sys.stderr.write("   "+key+": "+str(value))

    header="ref,pos,n,"
    if args.strand:
        header+=",".join([base+strand for base in bases for strand in strands])
    else:
        header+=','.join(bases)
    print(header)
    for column in countBasesInFile(args.bamFile,args.region):
        if args.strand:
            counts=[column[strand][base] for base in bases for strand in strands]
        else: 
            counts=[column['+'][base]+column['-'][base] for base in bases]
        out=','.join([str(x) for x in counts])
        print("%s,%d,%d,%s" % (column["ref"],column['pos'],column['n'],out))
            

    
if __name__ == '__main__':
    main(sys.argv[1:])
