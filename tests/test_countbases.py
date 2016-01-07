import pytest
from dnapy import countbases
import os
import stat
import pysam

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    with pytest.raises(ValueError):
        next(countbases.countBasesInFile(str(d)))
    #doesn't exist yet
    with pytest.raises(IOError):
        next(countbases.countBasesInFile(str(p)))
    #incorrectly formatted file
    p.write("test")
    with pytest.raises(ValueError):
        next(countbases.countBasesInFile(str(p)))
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(ValueError):
        next(countbases.countBasesInFile(str(p)))



def test_goodFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.bam')
    header = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'ref'}] }
    outFile=pysam.AlignmentFile(str(p),"wb",header=header)
    a = pysam.AlignedSegment()
    a.query_name = "read1"
    a.query_sequence="AAAAATTTTT"
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((0,10), )
    a.query_qualities = "(((((((((("
    outFile.write(a)
    outFile.close()
    pysam.index(str(p))
    count=0
    for col in countbases.countBasesInFile(str(p)):
        if count<5:
            assert col['A']==1
            assert col['G']+col['T']+col['C']==0
            assert col['n']==1
        else:
            assert col['T']==1
            assert col['G']+col['A']+col['C']==0
            assert col['n']==1
        assert col['pos']==count+32
        count+=1

    outFile=pysam.AlignmentFile(str(p),"wb",header=header)
    a = pysam.AlignedSegment()
    a.query_name = "read1"
    a.query_sequence="AAAAACCCCC"
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((0,10), )
    a.query_qualities = "(((((((((("
    outFile.write(a)
    a = pysam.AlignedSegment()
    a.query_name = "read2"
    a.query_sequence="AAAAATTTTT"
    a.reference_id = 0
    a.reference_start = 32
    a.mapping_quality = 20
    a.cigar = ((0,10), )
    a.query_qualities = "(((((((((("
    outFile.write(a)
    outFile.close()
    pysam.index(str(p))
    count=0
    for col in countbases.countBasesInFile(str(p)):
        if count<5:
            assert col['A']==2
            assert col['G']+col['T']+col['C']==0
            assert col['n']==2
        else:
            print(col)
            assert col['T']==1
            assert col['C']==1
            assert col['G']+col['A']==0
            assert col['n']==2
        assert col['pos']==count+32
        count+=1


