import os
import argparse
import gzip
import sys

def check_file(targetFile):
    if not os.path.isfile(targetFile):
        raise argparse.ArgumentTypeError(targetFile+' is not a file')
    if os.access(targetFile, os.R_OK):
        return targetFile
    else:
        raise argparse.ArgumentTypeError(targetFile+' is not readable')



def writeFastqRead(readFile,read):
    out="@%s\n%s\n+%s\n%s\n" %(read[0],read[1],read[0],read[2])
    out=out.encode('utf-8')
    readFile.write(out)


def closeFiles(openFiles):
    for openFile in openFiles:
        openFile.close()

def openNormalOrGz(gzFile,mode='r'):
    try:
        if gzFile[-2:]=='gz' or gzFile[-4]=='gzip':
            fastq=gzip.open(gzFile, mode)
        else:
            fastq=open(gzFile,mode)
    except IOError:
        sys.stderr.write("Problem opening file:"+gzFile+"\n")
        raise
    return fastq
