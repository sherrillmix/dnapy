import pytest
from dnapy import removereads
import os
import stat
import pysam
import subprocess

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    with pytest.raises(IOError):
        next(removereads.filterFastqIter([str(p)],set()))
    #doesn't exist yet
    with pytest.raises(IOError):
        removereads.filterFastqIter([str(p)],set())
    #incorrectly formatted file
    #doesn't exist yet
    p.write("@test\nAAAA\n+\n1") #qual and seq different
    with pytest.raises(ValueError):
        next(removereads.filterFastqIter([str(p)],set()))
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(IOError):
        removereads.filterFastqIter([str(p)],set())


def test_main(capsys,tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.fastq')
    p.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n((\n@seq3\nT\n+\n(\n")
    p2 = d.join('test2.fastq')
    p2.write("@seq1z\nTTT\n+\n(((\n@seq2z\nTT\n+\n((\n@seq3\nT\n+\n(\n")
    f = d.join('test.filter')
    f.write("seq2\nseq3")
    o = d.join('test.out')
    o2 = d.join('test2.out')
    removereads.main([str(p),'-f',str(f),'-o',str(o),'-d1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['.','Good reads: 1 Bad reads: 2']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in o.readlines()],['@seq1','AAA','+seq1','(((']):
        assert ii==jj
    removereads.main([str(p),str(p2),'-f',str(f),'-o',str(o)+","+str(o2),'-d1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['.','Good reads: 1 Bad reads: 2']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in o.readlines()],['@seq1','AAA','+seq1','(((']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in o2.readlines()],['@seq1z','TTT','+seq1z','(((']):
        assert ii==jj



