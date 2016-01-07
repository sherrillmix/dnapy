import argparse
import pysam
import collections
import os
import sys


MAX_DEPTH=1e9

def countBasesInFile(inputFile,region=None):
    samfile = pysam.AlignmentFile(inputFile, "rb" )
    print(samfile.pileup(region=region))
    for pileupcolumn in samfile.pileup(region=region,max_depth=MAX_DEPTH):
        #print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
        counts=collections.defaultdict(lambda:0)
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                thisBase=pileupread.alignment.query_sequence[pileupread.query_position]
                counts[thisBase]+=1
        yield {"ref": pileupcolumn.reference_name, "pos": pileupcolumn.pos, "n": pileupcolumn.n, "A": counts['A'], "C": counts['C'], "G": counts['G'], "T":counts['T']}
    samfile.close()



def main(argv=None):
    if argv is None:
        argv=sys.argv[1:]

    parser = argparse.ArgumentParser(description="A program to count the number of bases at each position in a region. The command generates a csv file (default: output.csv, use `-o` or `--output` to change) with 3 columns; start of region, end of region, , C, G, T.")
    parser.add_argument("-v","--verbose", help="increase output verbosity to stderr", action="store_true")
    parser.add_argument("-b","--bamFile", help="a bam file containing the alignment",type=check_file,required=True)
    parser.add_argument("-r","--region", help="the region to count in",type=string,default=None)
    args=parser.parse_args(argv)
        
    if args.verbose:
        sys.stderr.write("Arguments: ")
        for key, value in vars(args).items():
            sys.stderr.write("   "+key+": "+str(value))

        for column in countBasesInFile(args.bamFile,args.region):
            print ("%s,%d,%d,%d,%d,%d,%d" % (column["ref"],column['pos'],column['n'],column["A"],column['C'],column['G'],column["A"]))
            

    
if __name__ == '__main__':
    main(sys.argv[1:])
