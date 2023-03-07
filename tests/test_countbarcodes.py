import pytest
from dnapy import countbarcodes
from dnapy import helper
import os
import stat
import argparse

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p1 = d.join('test.txt')
    p2 = d.join('test1.txt')
    p3 = d.join('test2.txt')
    p2.write("@test\nAAAA\n+\n1111")
    p3.write("AAAA\nBBBB\nCCCC")
    #not a file
    with pytest.raises(IOError):
        next(countbarcodes.barcodeFastqIter([str(d)],1,4,str(p3)))
    #doesn't exist yet
    with pytest.raises(IOError):
        next(countbarcodes.barcodeFastqIter([str(p1)],1,4,str(p3)))
    #incorrectly formatted file
    p1.write("@test\nAAAA\n+\n1") #qual and seq different
    with pytest.raises(ValueError):
        next(countbarcodes.barcodeFastqIter([str(p1)],1,4,str(p3)))
    #make unreadable
    os.chmod(str(p1),os.stat(str(p1)).st_mode & ~stat.S_IREAD)
    with pytest.raises(IOError):
        next(countbarcodes.barcodeFastqIter([str(p1)],1,4,str(p3)))


def test_main(capsys,tmpdir):
    with pytest.raises(SystemExit):
        countbarcodes.main()
    out, err=capsys.readouterr()
    assert 'usage' in err
    with pytest.raises(SystemExit):
        countbarcodes.main(['-h'])
    out, err=capsys.readouterr()
    assert 'usage' in out

    d = tmpdir.mkdir('dir')
    p1 = d.join('test.fastq')
    p1.write("@seq1\nAAAAAAYAA\n+\n(((((((((\n@seq2\nTTTTTTZ\n+\n(AAAAAA\n@seq2\nTTTTTTG\n+\n(AAAAAA\n@seq3\nT\n+\n(\n")
    b = d.join('test.bar')
    b.write('AAAAAA\nTTTTTT\nZZZZZZ')
    bBad = d.join('test.bar2')
    o = d.join('test.out')
    with pytest.raises(SystemExit):
        #missing barcode file
        countbarcodes.main(['-w',str(bBad)])
    with pytest.raises(SystemExit):
        #missing fastq
        countbarcodes.main(['-w',str(b),'-s','1','-e','7'])

    out, err=capsys.readouterr()
    countbarcodes.main([str(p1),'-w',str(b),'-s','1','-e','6','-d','1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['...','Matching barcodes found: 3 Unassigned barcodes: 1']):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),['AAAAAA,1','TTTTTT,2']):
        assert ii==jj
    #Check no barcode list outputs all and reads shorter than barcode (currently output)
    countbarcodes.main([str(p1),'-s','1','-e','6','-d','1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['....','Matching barcodes found: 4 Unassigned barcodes: 0']):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),['AAAAAA,1','TTTTTT,2','T,1']):
        assert ii==jj
    #test no dots, no matches (due to short barcode def not matching input barcodes)
    countbarcodes.main([str(p1),'-s','1','-e','4'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),[]):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),[]):
        assert ii==jj
    #test output zeros
    countbarcodes.main([str(p1),'-w',str(b),'-s','1','-e','6','-d','1','-a'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['...','Matching barcodes found: 3 Unassigned barcodes: 1']):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),['AAAAAA,1','TTTTTT,2','ZZZZZZ,0']):
        assert ii==jj
    #test no whitelist and all option. currently no error, should there be?
    countbarcodes.main([str(p1),'-s','1','-e','6','-d','1','-a'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['....','Matching barcodes found: 4 Unassigned barcodes: 0']):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),['AAAAAA,1','TTTTTT,2','T,1']):
        assert ii==jj
    #test start end moves barcodes
    b.write('Y\nG\nZZZ') #currently no error if barcode doesn't match length, should there be?
    countbarcodes.main([str(p1),'-w',str(b),'-s','7','-e','7','-d','1'])
    out, err=capsys.readouterr()
    for ii,jj in zip(err.split('\n'),['..','Matching barcodes found: 2 Unassigned barcodes: 2']):
        assert ii==jj
    for ii,jj in zip(out.split('\n'),['Y,1','G,1']):
        assert ii==jj



