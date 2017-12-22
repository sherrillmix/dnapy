import os
import argparse
import gzip
import sys
import stat

def readSimpleCsv(csvFile):
    lines=[line.strip() for line in openNormalOrGz(csvFile)]
    splits=[[entry.strip('" \'') for entry in line.split(',')] for line in lines if line]
    n=[len(x) for x in splits]
    if any([x!=n[0] for x in n]):
        raise ValueError('All rows do not have same numbers of entries in '+csvFile)
    return splits

def readSimpleFastq(fileHandle):
    while True:
        try:
            seq=[next(fileHandle).strip() for _ in range(4)]
            del(seq[2])
            seq[0]=seq[0][1:]
            yield seq
        except StopIteration:
            raise

def checkFile(targetFile,orPipe=True):
    if not os.path.exists(targetFile):
            raise argparse.ArgumentTypeError(targetFile+' does not exist')
    if not os.path.isfile(targetFile):
        if not orPipe or not stat.S_ISFIFO(os.stat(targetFile).st_mode):
            raise argparse.ArgumentTypeError(targetFile+' is not a file')
    if os.access(targetFile, os.R_OK):
        return targetFile
    else:
        raise argparse.ArgumentTypeError(targetFile+' is not readable')

def checkDir(targetFile,create=True):
    if create and not os.path.exists(targetFile):
        os.makedirs(targetFile)
    if not os.path.isdir(targetFile):
        raise argparse.ArgumentTypeError(targetFile+' is not a directory')
    if os.access(targetFile, os.W_OK):
        return targetFile
    else:
        raise argparse.ArgumentTypeError(targetFile+' is not writeable')

def writeFastqRead(readFile,read):
    out="@%s\n%s\n+%s\n%s\n" %(read[0],read[1],read[0],read[2])
    #out=out.encode('utf-8')
    readFile.write(out)

def closeFiles(openFiles):
    try:
        iterator = openFiles.items()
    except AttributeError:
        iterator = enumerate(openFiles)
    for _, openFile in iterator:
        openFile.close()

def openNormalOrGz(gzFile,mode='r'):
    if not any([True for ii in mode if ii=='t']):mode+="t"
    try:
        if gzFile[-2:]=='gz' or gzFile[-4]=='gzip':
            fastq=gzip.open(gzFile, mode)
        else:
            fastq=open(gzFile,mode)
    except IOError:
        sys.stderr.write("Problem opening file:"+gzFile+"\n")
        raise
    return fastq

def checkPositiveInt(value,maxVal=None):
    ivalue=int(value)
    if ivalue<=0:
        raise argparse.ArgumentTypeError('%s is not a postive integer' % value)
    if maxVal is not None:
        if ivalue > maxVal: 
            raise argparse.ArgumentTypeError('%s is not less than %d' % (value,maxVal))
    return ivalue

