import argparse
import pysam
import sys
from dnapy import helper



def getAlignsInFile(inputFile,region=None,minQuality=0,endspan=0):
    with pysam.AlignmentFile(inputFile, "rb" ) as samfile:
        for read in samfile.fetch(region=region):
            if read.cigartuples is None or read.mapping_quality<minQuality:
                continue
            rPos=read.reference_start
            qPos=0
            tPos=read.pos
            qPos=0
            seq=''
            insertions=[]
            cigarOps=[[operation,length] for operation, length in read.cigartuples]
            if endspan>0:
                for ii in range(len(cigarOps)-1,-1,-1):
                    if cigarOps[ii][0] in [0,7,8]:
                        if cigarOps[ii][1]>endspan:
                            break
                        cigarOps[ii][0]=4
                for ii in range(0,len(cigarOps)):
                    if cigarOps[ii][0] in [0,7,8]:
                        if cigarOps[ii][1]>endspan:
                            break
                        cigarOps[ii][0]=4
            for operation,length in cigarOps:
                if operation in [0,7,8]:
                    seq+=read.query_sequence[qPos:(qPos+length)]
                if operation in [2,3]:
                    seq+='-'*length
                if operation in [1]:
                    insertions.append([tPos,read.query_sequence[qPos:(qPos+length)]])
                if operation in [0,2,3,7,8]:
                    tPos+=length
                if operation in [0,1,4,7,8]:
                    qPos+=length
                #5,6 just ignore
                if operation < 0 or operation > 8:
                    raise ValueError('Unknown operation in read "'+read.query_name+'" cigar: '+read.cigarstring)
            # looks like sequence is already reverse complimented if strand is -
            yield {"name": read.query_name, "start": read.pos, "seq": seq,'strand': '-' if read.flag&16==16 else '+','insertions': insertions}

def padRead(seq,start,end,readInserts=None,refInserts=None):
    complete='-'*start+seq+'-'*(end-len(seq)-start)
    if len(seq)+start>end: raise IndexError('Start position plus read length greater than end')
    if refInserts is None: return complete
    inserts={pos: ['-']*length for (pos,length) in refInserts.iteritems()}
    if readInserts is not None:
        for pos,seq in readInserts:
            if pos not in inserts: raise IndexError('Read insert not in total ref inserts')
            if len(seq)>len(inserts[pos]): raise IndexError('Read insert longer than ref insert')
            inserts[pos][0:len(seq)]=seq
    sortInserts=[(pos,''.join(seqArray)) for (pos,seqArray) in inserts.iteritems()]
    sortInserts.sort(key=lambda x: x[0])
    splitPos=[xx[0] for xx in sortInserts]
    seqSplit=[complete[start:end] for start,end in zip([0]+splitPos,splitPos+[None])]
    return ''.join([val for pair in zip(['']+[xx[1] for xx in sortInserts],seqSplit) for val in pair])

def getRefFromFasta(fasta,region):
    ref=helper.readSimpleFasta(fasta)
    if len(ref)==0:
        raise ImportError('Fasta file is empty')
    if len(ref)>1:
        if region is None:
            raise ImportError('Reference fasta file contains more than one reference and region unspecified')
    if region is None:
        ref=ref[0]
        region=ref[0].split(' ',1)[0]
        ref=ref[1]
    else:
        target=region.split(':',1)[0]
        refIds=[ii for ii in range(len(ref)) if ref[ii][0].split(' ',1)[0]==target]
        if len(refIds)==0:
            raise ImportError('Reference fasta file does not contain reference matching region '+region)
        if len(refIds)>1:
            raise ImportError('Reference fasta file contains multiple references matching region '+region)
        ref=ref[refIds[0]][1]
    return [ref,region]



def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to convert a bam file into an aligned fasta file. The command generates fasta formatted output (two lines for each sequence: a name line prepended by > and a line containing the aligned sequence) to standard out.")
    parser.add_argument('bamFile', help='a bam file containing the alignment',type=helper.checkFile)
    parser.add_argument("-s","--refseq", help="fasta file giving the reference sequence of interest",type=helper.checkFile,required=True)
    parser.add_argument("-q","--minQuality", help="don't count alignments with a mapping quality less than this", type=int,default=0)
    parser.add_argument("-v","--verbose", help="increase output verbosity to stderr", action="store_true")
    parser.add_argument("-r","--region", help="the region to pull reads from",default=None)
    parser.add_argument("-e","--endspan", help="ignore spans of matches at the start or end of a read less than this cutoff",default=0,type=int)
    args=parser.parse_args(argv)
 
        
    if args.verbose:
        sys.stderr.write("Arguments: \n")
        for key, value in vars(args).items():
            sys.stderr.write("   "+key+": "+str(value)+'\n')

    with helper.openNormalOrGz(args.refseq) as fasta:
        ref,args.region=getRefFromFasta(fasta,args.region)

    nRead=0
    aligns=[read for read in getAlignsInFile(args.bamFile,args.region,args.minQuality,args.endspan)]
    inserts=[[insertion[0],len(insertion[1])] for align in aligns if len(align['insertions'])>0 for insertion in align['insertions'] ]
    maxInserts={}
    for pos,length in inserts:
        maxInserts[pos]=max(length,maxInserts[pos]) if pos in maxInserts else length
    nChar=len(ref)
    print('>'+args.region)
    refPad=padRead(ref,0,nChar,[],maxInserts)
    print(refPad)
    for align in aligns:
        padded=padRead(align['seq'],align['start'],nChar,align['insertions'],maxInserts)
        print('>'+align['name'])
        print(padded)
    
if __name__ == '__main__':
    main()
