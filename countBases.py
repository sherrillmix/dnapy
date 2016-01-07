import pysam
from collections import defaultdict


def check_file(targetFile):
    if not os.path.isfile(targetFile):
        raise argparse.ArgumentTypeError(targetFile+' is not a file')
    if os.access(targetFile, os.R_OK):
        return targetFile
    else:
        raise argparse.ArgumentTypeError(targetFile+' is not readable')


def main(argv=None):
    if argv is None:
        argv=sys.argv[1:]

    parser = argparse.ArgumentParser(description="A program to count the number of bases at each position in a region. The command generates a csv file (default: output.csv, use `-o` or `--output` to change) with 3 columns; start of region, end of region, , C, G, T.")
    parser.add_argument("-v","--verbose", help="increase output verbosity to stderr", action="store_true")
    parser.add_argument("-b","--bamFile", help="a bam file containing the alignment",type=check_file,required=True)
    parser.add_argument("-r","--region", help="the region to count in",type=string,default=None)
    args=parser.parse_args(argv)
        
    if args.verbose:
        print("Arguments: ")
        for key, value in vars(args).items():
            print("   "+key+": "+str(value))

samfile = pysam.AlignmentFile("S9R1_S13vH3N2_HA_Colorado_edited.bam", "rb" )
for pileupcolumn in samfile.pileup("H3N2_HA_Colorado_edited",max_depth=1e9):
    #print ("coverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
    counts=defaultdict(lambda:0)
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            thisBase=pileupread.alignment.query_sequence[pileupread.query_position]
            counts[thisBase]+=1
    print ("%d,%d,%d,%d,%d,%d" % (pileupcolumn.pos, pileupcolumn.n, counts['A'], counts['C'], counts['G'], counts['T']))

samfile.close()
