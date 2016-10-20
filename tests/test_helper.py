import pytest
import os
import stat
from dnapy import helper
import argparse


def test_check_file(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    with pytest.raises(argparse.ArgumentTypeError):
        helper.check_file(str(d))
    #doesn't exist yet
    with pytest.raises(argparse.ArgumentTypeError):
        helper.check_file(str(p))
    p.write("test")
    assert helper.check_file(str(p))==str(p)
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(argparse.ArgumentTypeError):
        helper.check_file(str(p))

def test_check_dir(tmpdir):
    d = tmpdir.mkdir('dir')
    assert helper.check_dir(str(d))==str(d)
    p = d.join('test.txt')
    with pytest.raises(argparse.ArgumentTypeError):
        helper.check_dir(str(p))
    p.write("test")
    with pytest.raises(argparse.ArgumentTypeError):
        helper.check_dir(str(p))
    #make unwriteable
    os.chmod(str(d),os.stat(str(d)).st_mode & ~stat.S_IWRITE)
    with pytest.raises(argparse.ArgumentTypeError):
        helper.check_dir(str(d))

def test_openGzOrNormal(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    gzFile = d.join('test.gz')
    gz=helper.openNormalOrGz(str(gzFile),'w')
    helper.writeFastqRead(gz,['1','22','333'])
    helper.writeFastqRead(gz,['55555','666666','7777777'])
    helper.closeFiles([gz])
    gz=helper.openNormalOrGz(str(gzFile))
    pred=['@1','22','+1','333','@55555','666666','+55555','7777777']
    for x,y in zip(gz,[x+'\n' for x in pred]):
        assert x==y


