import pytest
from dnapy import splitbarcodes
from dnapy import helper
import os
import stat
import argparse

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
    o = d.join('test_1.fastq.gz')
    o2 = d.join('test_2.fastq.gz')
    #duplicate barcode
    b.write("test,AAA\ntest2,AAA\ntest3,AAT")
    with pytest.raises(argparse.ArgumentTypeError):
        splitbarcodes.main([str(p1),'-i',str(p1),'-b',str(b),'-o',str(d),'-d1'])
    #duplicate sample name
    b.write("test,AAA\ntest,AAC\ntest3,AAT")
    with pytest.raises(argparse.ArgumentTypeError):
        splitbarcodes.main([str(p1),'-i',str(p1),'-b',str(b),'-o',str(d),'-d1'])
    #sample named __UNASSIGNED__R
    b.write("__UNASSIGNED__R,AAA")
    with pytest.raises(argparse.ArgumentTypeError):
        splitbarcodes.main([str(p1),'-i',str(p1),'-b',str(b),'-o',str(d),'-d1','-u'])
    b.write("__UNASSIGNED__I,AAA")
    with pytest.raises(argparse.ArgumentTypeError):
        splitbarcodes.main([str(p1),'-i',str(p1),'-b',str(b),'-o',str(d),'-d1','-u'])
    b.write("test,AAA")
    #two index files, one barcode
    with pytest.raises(argparse.ArgumentTypeError):
        splitbarcodes.main([str(p1),'-i',str(p1),str(p1),'-b',str(b),'-o',str(d),'-d1'])

    #one read file, one index file
    splitbarcodes.main([str(p1),'-i',str(p1),'-b',str(b),'-o',str(d),'-d1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['.','Reads assigned to barcode: 1 Unassigned reads: 2']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o)).readlines()],['@seq1','AAA','+seq1','(((']):
        assert ii==jj
    
    #two read files, one index file
    splitbarcodes.main([str(p1),str(p2),'-i',str(p1),'-b',str(b),'-o',str(d),'-d1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['.','Reads assigned to barcode: 1 Unassigned reads: 2']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o)).readlines()],['@seq1','AAA','+seq1','(((']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o2)).readlines()],['@seq1z','TTT','+seq1z','(((']):
        assert ii==jj

    #two read, two index files. test unassigned
    b.write("test,T,T\ntest2,A,T")
    i1 = d.join('test3.fastq')
    i1.write("@seq1\nT\n+\n(\n@seq2\nA\n+\n(\n@seq3\nT\n+\n(\n")
    i2 = d.join('test4.fastq')
    i2.write("@seq1z\nT\n+\n(\n@seq2z\nT\n+\n(\n@seq3\nTT\n+\n(!\n")
    splitbarcodes.main([str(p1),str(p2),'-i',str(i1),str(i2),'-b',str(b),'-o',str(d),'-d1','-u'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['..','Reads assigned to barcode: 2 Unassigned reads: 1']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o)).readlines()],['@seq1','AAA','+seq1','(((','@seq2','TT','+seq2','(A']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(o2)).readlines()],['@seq1z','TTT','+seq1z','(((','@seq2z','TT','+seq2z','((']):
        assert ii==jj
    #test the unassigned
    ur1=d.join('__UNASSIGNED__R1.fastq.gz')
    ur2=d.join('__UNASSIGNED__R2.fastq.gz')
    ui1=d.join('__UNASSIGNED__I1.fastq.gz')
    ui2=d.join('__UNASSIGNED__I2.fastq.gz')
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(ur1)).readlines()],['@seq3','T','+seq3','(']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(ur2)).readlines()],['@seq3','TAA','+seq3','(!!']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(ui1)).readlines()],['@seq3','T','+seq3','(']):
        assert ii==jj
    for ii,jj in zip([x.rstrip('\n') for x in helper.openNormalOrGz(str(ui2)).readlines()],['@seq3','TT','+seq3','(!']):
        assert ii==jj



