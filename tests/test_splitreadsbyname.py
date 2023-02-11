import pytest
from dnapy import splitreadsbyname
from dnapy import helper
import os
import stat
import argparse
import Bio

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    #not a file
    with pytest.raises(IOError):
        next(splitreadsbyname.splitFastqIter(str(d),set()))
    #doesn't exist yet
    with pytest.raises(IOError):
        splitreadsbyname.splitFastqIter(str(p),set())
    #incorrectly formatted file
    p.write("@test\nAAAA\n+\n1") #qual and seq different
    with pytest.raises(ValueError):
        next(splitreadsbyname.splitFastqIter(str(p),set()))
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(IOError):
        splitreadsbyname.splitFastqIter(str(p),set())



def test_main(capsys,tmpdir):
    with pytest.raises(SystemExit):
        splitreadsbyname.main()
    out, err=capsys.readouterr()
    assert 'usage' in err
    with pytest.raises(SystemExit):
        splitreadsbyname.main(['-h'])
    out, err=capsys.readouterr()
    assert 'usage' in out

    d = tmpdir.mkdir('dir')
    p = d.join('test.fastq')
    p.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n((\n@seq3\nT\n+\n(\n@seq4\nC\n+\n3\n@seq5\nA\n+\n1\n@seq6 OTHER\nZZ\n+\n12\n@seq6\nYY\n+\nZZ")
    f = d.join('test.filter')
    f.write("seq2,a\nseq3,a\nseq99,z\nseq6,bb")
    o = d.mkdir('test.out')
    with pytest.raises(SystemExit): splitreadsbyname.main([str(p),str(p),'-s',str(f),'-o',str(o),'-d1'])
    with pytest.raises(SystemExit):
        splitreadsbyname.main(['-s',str(f),'-o',str(o),'-d1'])
    with pytest.raises(SystemExit):
        splitreadsbyname.main([str(p),'-s','-o',str(o),'-d1'])
    with pytest.raises(SystemExit):
        splitreadsbyname.main([str(p),'-o',str(o),'-d1'])
    #clear 
    out, err=capsys.readouterr()
    splitreadsbyname.main([str(p),'-s',str(f),'-o',str(o),'-d1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['....','Reads assigned to files: 4 Unassigned reads: 3']): assert ii==jj
    for ii,jj,kk in zip(Bio.SeqIO.QualityIO.FastqGeneralIterator(o.join('a')),['seq2','seq3'],['TT','T']):
        assert ii[0]==jj
        assert ii[1]==kk
    #duplicate read
    f.write("seq2,a\nseq3,a\nseq99,z\nseq2,c") 
    with pytest.raises(argparse.ArgumentTypeError):
        splitreadsbyname.main([str(p),'-s',str(f),'-o',str(o),'-d1'])
    #make sure doesn't append
    f.write("seq2,a\nseq3,a\nseq99,z")
    p.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n((\n@seq3\nT\n+\n(\n@seq4\nC\n+\n3\n@seq5\nA\n+\n1")
    splitreadsbyname.main([str(p),'-s',str(f),'-o',str(o),'-d1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['..','Reads assigned to files: 2 Unassigned reads: 3']): assert ii==jj
    for ii,jj,kk in zip(Bio.SeqIO.QualityIO.FastqGeneralIterator(o.join('a')),['seq2','seq3'],['TT','T']):
        assert ii[0]==jj
        assert ii[1]==kk
    #make sure does append
    splitreadsbyname.main([str(p),'-s',str(f),'-o',str(o),'-d1','-a','-u'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['..','Reads assigned to files: 2 Unassigned reads: 3']): assert ii==jj
    for ii,jj,kk in zip(Bio.SeqIO.QualityIO.FastqGeneralIterator(o.join('a')),['seq2','seq3','seq2','seq3'],['TT','T','TT','T']):
        assert ii[0]==jj
        assert ii[1]==kk
    for ii,jj,kk in zip(Bio.SeqIO.QualityIO.FastqGeneralIterator(helper.openNormalOrGz(o.join('__UNASSIGNED__.fastq.gz'))),['seq1','seq4','seq5'],['AAA','C','A']):
        assert ii[0]==jj
        assert ii[1]==kk
 
