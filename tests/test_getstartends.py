import pytest
from dnapy import getstartends
import os
import stat
import pysam
import subprocess

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    with pytest.raises(ValueError):
        next(getstartends.getStartsInFile(str(d)))
    #doesn't exist yet
    with pytest.raises(IOError):
        next(getstartends.getStartsInFile(str(p)))
    #incorrectly formatted file
    p.write("test")
    with pytest.raises(ValueError):
        next(getstartends.getStartsInFile(str(p)))
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(ValueError):
        next(getstartends.getStartsInFile(str(p)))

#something about capsys messes up pysam.index so create ahead of time for testing main
@pytest.fixture(scope='session')
def bamFile(tmpdir_factory):
    header = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'ref'}] }
    p=tmpdir_factory.mktemp('test').join('test.bam')
    outFile=pysam.AlignmentFile(str(p),"wb",header=header)
    a = pysam.AlignedSegment()
    a.query_name = "read3"
    a.query_sequence="GGGGAAAAAT"
    a.reference_start = 28
    a.reference_id = 0
    a.mapping_quality = 20
    a.cigar = ((0,10), )
    #a.query_qualities = pysam.qualitystring_to_array("((((((((((")
    a.flag=16
    outFile.write(a)
    a.query_name = "read2"
    a.reference_start = 32
    a.query_sequence="AAAAATTTTT"
    a.flag=0
    outFile.write(a)
    a.query_name = "read1"
    a.query_sequence="TTAAAAACCCCCGGC"
    #a.query_qualities = pysam.qualitystring_to_array("(((((((((((((")
    a.cigar = ((5,5),(4,2),(0,10), (2,2),(0,1),(1,1),(0,1))
    outFile.write(a)
    outFile.close()
    pysam.index(str(p))
    return(p)

def test_goodFiles(tmpdir,bamFile):
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
    #a.query_qualities = pysam.qualitystring_to_array("((((((((((")
    outFile.write(a)
    outFile.close()
    pysam.index(str(p))
    count=0
    out=next(getstartends.getStartsInFile(str(p)))
    assert out['start']==33
    assert out['end']==42
    assert out['strand']=='+'
    assert out['ref']=='ref'
    for read,start,strand,end in zip(getstartends.getStartsInFile(str(bamFile),maxGaps=10),[29,33,33],['-','+','+'],[38,42,46]):
        assert read['start']==start
        assert read['strand']==strand
        assert read['end']==end
    for read,start,strand,end in zip(getstartends.getStartsInFile(str(bamFile)),[29,33],['-','+'],[38,42]):
        assert read['start']==start
        assert read['strand']==strand
        assert read['end']==end


def test_main(capsys,tmpdir,bamFile):
    with pytest.raises(SystemExit):
        getstartends.main()
    out, err=capsys.readouterr()
    assert 'usage' in err
    with pytest.raises(SystemExit):
        getstartends.main(['-h'])
    out, err=capsys.readouterr()
    assert 'usage' in out

    getstartends.main(['-v',str(bamFile)])
    out, err=capsys.readouterr()
    assert 'Arguments' in err
    compare=getstartends.getStartsInFile(str(bamFile))

    count=0
    for ii,jj in zip(out.split('\n')[1:],compare): 
        ii=ii.split(',')
        count+=1
        assert ii[0]==jj['ref']
        assert int(ii[1])==jj['start']
        assert int(ii[2])==jj['end']
        assert ii[3]==jj['strand']
    assert count==2

    compare=getstartends.getStartsInFile(str(bamFile),maxGaps=10)
    getstartends.main([str(bamFile),'-g 10'])
    out, err=capsys.readouterr()

    count=0
    for ii,jj in zip(out.split('\n')[1:],compare): 
        count+=1
        ii=ii.split(',')
        assert ii[0]==jj['ref']
        assert int(ii[1])==jj['start']
        assert int(ii[2])==jj['end']
        assert ii[3]==jj['strand']
    assert count==3

    getstartends.main([str(bamFile),'-g 10', '-n'])
    noHead, err=capsys.readouterr()
    for ii,jj in zip(out.split('\n')[1:],noHead.split('\n')): 
        assert ii==jj

    getstartends.main([str(bamFile),'-g 10', '-n','-c'])
    regCol, err=capsys.readouterr()
    for ii,jj in zip(out.split('\n')[1:],regCol.split('\n')):
        if ii:
            ii="None,"+ii
        assert ii==jj

    getstartends.main([str(bamFile), '-g 10','-n','-c','-r', 'ref:1-100'])
    regCol, err=capsys.readouterr()
    for ii,jj in zip(out.split('\n')[1:],regCol.split('\n')):
        if ii:
            ii="ref:1-100,"+ii
        assert ii==jj


    d = tmpdir.mkdir('dir')
    p = d.join('test.regions')
    p.write("ref:1-10\nref:1-100")
    getstartends.main([str(bamFile), '-g 10','-n','-c','-f', str(p)])
    fileOut, err=capsys.readouterr()
    for ii,jj in zip(regCol.split('\n'),fileOut.split('\n')):
        assert ii==jj


def test_commandline(capsys,bamFile):
    getstartends.main([str(bamFile)])
    out, err=capsys.readouterr()
    out2 = subprocess.check_output("getstartends "+str(bamFile), shell=True)
    out3 = subprocess.check_output("python -m dnapy.getstartends "+str(bamFile), shell=True)
    for ii,jj,kk in zip(out.split('\n'),out2.decode().split('\n'),out3.decode().split('\n')): 
        print(ii)
        print(jj)
        assert ii==jj==kk
