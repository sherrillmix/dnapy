import pytest
from dnapy import splitbarcodes
from dnapy import helper
import os
import stat

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p1 = d.join('test.txt')
    p2 = d.join('test1.txt')
    p2.write("@test\nAAAA\n+\n1111")
    #not a file
    with pytest.raises(IOError):
        next(splitbarcodes.barcodeFastqIter([str(d)],[str(p2)],set()))
    with pytest.raises(IOError):
        next(splitbarcodes.barcodeFastqIter([str(p2)],[str(d)],set()))
    #doesn't exist yet
    with pytest.raises(IOError):
        splitbarcodes.barcodeFastqIter([str(p1)],[str(p2)],set())
    with pytest.raises(IOError):
        splitbarcodes.barcodeFastqIter([str(p2)],[str(p1)],set())
    #incorrectly formatted file
    p1.write("@test\nAAAA\n+\n1") #qual and seq different
    with pytest.raises(ValueError):
        next(splitbarcodes.barcodeFastqIter([str(p1)],[str(p2)],set()))
    with pytest.raises(ValueError):
        next(splitbarcodes.barcodeFastqIter([str(p2)],[str(p1)],set()))
    #make unreadable
    os.chmod(str(p1),os.stat(str(p1)).st_mode & ~stat.S_IREAD)
    with pytest.raises(IOError):
        splitbarcodes.barcodeFastqIter([str(p1)],[str(p2)],set())
    with pytest.raises(IOError):
        splitbarcodes.barcodeFastqIter([str(p2)],[str(p1)],set())


def test_main(capsys,tmpdir):
    with pytest.raises(SystemExit):
        splitbarcodes.main()
    out, err=capsys.readouterr()
    assert 'usage' in err
    with pytest.raises(SystemExit):
        splitbarcodes.main(['-h'])
    out, err=capsys.readouterr()
    assert 'usage' in out

    d = tmpdir.mkdir('dir')
    p1 = d.join('test.fastq')
    p1.write("@seq1\nAAA\n+\n(((\n@seq2\nTT\n+\n(A\n@seq3\nT\n+\n(\n")
    p2 = d.join('test2.fastq')
    p2.write("@seq1z\nTTT\n+\n(((\n@seq2z\nTT\n+\n((\n@seq3\nT\n+\n(\n")
    b = d.join('test.filter')
    b.write("test,AAA")
    o = d.join('test_1.fastq.gz')
    o2 = d.join('test_2.fastq.gz')
    splitbarcodes.main([str(p1),'-i',str(p1),'-b',str(b),'-o',str(d),'-d1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['.','Good reads: 1 Bad reads: 2']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o)).readlines()],['@seq1','AAA','+seq1','(((']):
        assert ii==jj
    splitbarcodes.main([str(p1),str(p2),'-i',str(p1),'-b',str(b),'-o',str(d),'-d1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['.','Good reads: 1 Bad reads: 2']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o)).readlines()],['@seq1','AAA','+seq1','(((']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o2)).readlines()],['@seq1z','TTT','+seq1z','(((']):
        assert ii==jj
    b.write("test,T,T")
    i1 = d.join('test3.fastq')
    i1.write("@seq1\nT\n+\n(\n@seq2\nT\n+\n(\n@seq3\nT\n+\n(\n")
    i2 = d.join('test4.fastq')
    i2.write("@seq1z\nT\n+\n(\n@seq2z\nT\n+\n(\n@seq3\nTT\n+\n((\n")
    splitbarcodes.main([str(p1),str(p2),'-i',str(i1),str(i2),'-b',str(b),'-o',str(d),'-d1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['..','Good reads: 2 Bad reads: 1']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o)).readlines()],['@seq1','AAA','+seq1','(((','@seq2','TT','+seq2','(A']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o2)).readlines()],['@seq1z','TTT','+seq1z','(((','@seq2z','TT','+seq2z','((']):
        assert ii==jj

