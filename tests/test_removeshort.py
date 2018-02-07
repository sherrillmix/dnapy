import pytest
from dnapy import removeshort
from dnapy import helper
import os
import stat
import pysam
import subprocess

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    with pytest.raises(IOError):
        next(removeshort.shortFilterFastqIter(str(p)))
    #doesn't exist yet
    with pytest.raises(IOError):
        removeshort.shortFilterFastqIter(str(p))
    #incorrectly formatted file
    #doesn't exist yet
    p.write("@test\nAAAA\n+\n1") #qual and seq different
    with pytest.raises(ValueError):
        next(removeshort.shortFilterFastqIter(str(p)))
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(IOError):
        removeshort.shortFilterFastqIter(str(p))


def test_main(capsys,tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.fastq')
    p.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n((\n@seq3\nT\n+\n(\n")
    removeshort.main([str(p),'-l 2','-d 1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['..','Good reads: 2 Bad reads: 1']):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),['@seq1','AAA','+seq1','(((','@seq2','TT','+seq2','((']):
        assert ii==jj
    p.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n((\n@seq3\nT\n+\n(\n@seq4\nTTTTN\n+\n(((((")
    removeshort.main([str(p),'-l 2','-d 1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['...','Good reads: 3 Bad reads: 1']):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),['@seq1','AAA','+seq1','(((','@seq2','TT','+seq2','((','@seq4','TTTTN','+seq4','(((((']):
        assert ii==jj
    removeshort.main([str(p),'-l 2','-d 1','-n'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['..','Good reads: 2 Bad reads: 2']):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),['@seq1','AAA','+seq1','(((','@seq2','TT','+seq2','((']):
        assert ii==jj
    p.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n((\n@seq3\nTTT\n+\n(\n@seq4\nTTTTN\n+\n(((((")
    with pytest.raises(ValueError):
        removeshort.main([str(p),'-l 2','-d 1'])
    out, err=capsys.readouterr()
    p.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n((\n@seq3\nTTT\n+\n(\n@seq4\nTTTTN\n+\n(((((\n@seq5\nAA\n+\n(((")
    removeshort.main([str(p),'-l 2','-d 1','-p','-n'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['..','Good reads: 2 Bad reads: 3']):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),['@seq1','AAA','+seq1','(((','@seq2','TT','+seq2','((']):
        assert ii==jj
    badf=d.join('bad.fastq')
    removeshort.main([str(p),'-l 2','-d 1','-p','-n','-b',str(badf)])
    with open(str(badf)) as badh:
        for ii,jj in zip(badh.read().splitlines(),['@seq3','TTT','+seq3','(','@seq5','AA','+seq5','(((']):
            assert ii==jj


