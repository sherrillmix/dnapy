import pytest
from dnapy import helper
from dnapy import countkmers
import os
import stat
import collections

def test_splitIntoKmer():
    assert countkmers.splitIntoKmer("1234567",1)==['1','2','3','4','5','6','7']
    assert countkmers.splitIntoKmer("1234567",1,2)==['1','3','5','7']
    assert countkmers.splitIntoKmer("1234567",1,3)==['1','4','7']
    assert countkmers.splitIntoKmer("1234567",2,3)==['12','45']
    assert countkmers.splitIntoKmer("1234567",3)==['123','456']
    assert countkmers.splitIntoKmer("1234567",2)==['12','34','56']
    assert countkmers.splitIntoKmer("1234567",5)==['12345']
    assert countkmers.splitIntoKmer("1234567",5,100)==['12345']
    assert countkmers.splitIntoKmer("1234567",90)==[]
    
def test_generateAllKmers():
    assert countkmers.generateAllKmers(0)==[]
    assert countkmers.generateAllKmers(1)==['A','C','G','T']
    assert countkmers.generateAllKmers(1,['Z','y'])==['Z','y']
    assert countkmers.generateAllKmers(2,['Z','y'])==['ZZ','Zy','yZ','yy']
    assert len(countkmers.generateAllKmers(6,['A','B','C','D','E']))==pow(5,6)
    assert len(countkmers.generateAllKmers(10))==pow(4,10)
    assert len(set(countkmers.generateAllKmers(10)))==pow(4,10)

def test_countKmersInReads():
    out=collections.defaultdict(lambda:0)
    out['AT']=1
    assert countkmers.countKmersInReads([[1,'AT']],2)==out
    out['AT']=3
    assert countkmers.countKmersInReads([[1,'ATATAT']],2)==out
    out['CC']=2
    assert countkmers.countKmersInReads([[1,'ATATATCCCC']],2)==out
    assert countkmers.countKmersInReads([[1,'ATAT'],[1,'AT'],[2,'CCCC']],2)==out
    out=collections.defaultdict(lambda:0)
    out['ATAT']=1
    out['CCCC']=1
    assert countkmers.countKmersInReads([[1,'ATATCCCCZZZ']],4)==out


def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p1 = d.join('test.txt')
    p2 = d.join('test1.txt')
    p2.write("@test\nAAAA\n+\n1111")
    #not a file
    with pytest.raises(SystemExit):
        next(countkmers.main([str(d)]))
    with pytest.raises(SystemExit):
        next(countkmers.main([str(d),str(p2)]))
    with pytest.raises(SystemExit):
        next(countkmers.main([str(p2),str(d)]))
    #doesn't exist yet
    with pytest.raises(SystemExit):
        countkmers.main([str(p1)])
    with pytest.raises(SystemExit):
        countkmers.main([str(p1),str(p2)])
    with pytest.raises(SystemExit):
        countkmers.main([str(p2),str(p1)])
    #incorrectly formatted file
    p1.write("@test\nAAAA\n+\n1") #qual and seq different
    with pytest.raises(ValueError):
        next(countkmers.main([str(p1)]))
    with pytest.raises(ValueError):
        next(countkmers.main([str(p1),str(p2)]))
    with pytest.raises(ValueError):
        next(countkmers.main([str(p2),str(p1)]))
    #make unreadable
    os.chmod(str(p1),os.stat(str(p1)).st_mode & ~stat.S_IREAD)
    with pytest.raises(SystemExit):
        countkmers.main([str(p1)])
    with pytest.raises(SystemExit):
        countkmers.main([str(p2),str(p1)])
    with pytest.raises(SystemExit):
        countkmers.main([str(p1),str(p2)])


def test_main(capsys,tmpdir):
    with pytest.raises(SystemExit):
        countkmers.main()
    out, err=capsys.readouterr()
    assert 'usage' in err
    with pytest.raises(SystemExit):
        countkmers.main(['-h'])
    out, err=capsys.readouterr()
    assert 'usage' in out
    d = tmpdir.mkdir('dir')
    p1 = d.join('test1.fastq')
    p2 = d.join('test2.fastq')
    p3 = d.join('test3.fastq')
    p1.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n(A\n@seq3\nT\n+\n(\n")
    p2.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n(A\n@seq3\nT\n+\n(\n")
    p3.write("@seq1\nAAAA\n+\n((((\n@seq2\nCCC\n+\n(A(\n@seq3\nT\n+\n(\n")
    countkmers.main([str(p1),'-k2'])
    out, err=capsys.readouterr()
    for ii,jj in zip(out.split('\n'),['kmer,'+str(p1),'AA,1','TT,1']):
        assert ii==jj
    countkmers.main([str(p1),str(p2),'-k2'])
    out, err=capsys.readouterr()
    for ii,jj in zip(out.split('\n'),['kmer,'+str(p1)+','+str(p2),'AA,1,1','TT,1,1']):
        assert ii==jj
    countkmers.main([str(p1),str(p3),'-k2'])
    out, err=capsys.readouterr()
    for ii,jj in zip(out.split('\n'),['kmer,'+str(p1)+','+str(p3),'AA,1,2','CC,0,1','TT,1,0']):
        assert ii==jj
    countkmers.main([str(p1),str(p3),'-k3'])
    out, err=capsys.readouterr()
    for ii,jj in zip(out.split('\n'),['kmer,'+str(p1)+','+str(p3),'AAA,1,1','CCC,0,1']):
        assert ii==jj


