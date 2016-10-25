import pytest
import os
import stat
from dnapy import helper
import argparse


def test_readSimpleCsv(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    with pytest.raises(argparse.ArgumentTypeError):
        helper.checkFile(str(p))
    with helper.openNormalOrGz(str(p),'w') as f:
        f.write("1,'2',3   \n  \n  \na,\"bb\",ccc\n  2,3,4  ")
    assert(helper.readSimpleCsv(str(p))==[['1','2','3'],['a','bb','ccc'],["2","3","4"]])
    with helper.openNormalOrGz(str(p),'w') as f:
        f.write("1,2,3\n\n  \na,bb,ccc,d\n  2,3,4  ")
    with pytest.raises(ValueError):
        helper.readSimpleCsv(str(p))

def test_checkFile(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    with pytest.raises(argparse.ArgumentTypeError):
        helper.checkFile(str(d))
    #doesn't exist yet
    with pytest.raises(argparse.ArgumentTypeError):
        helper.checkFile(str(p))
    p.write("test")
    assert helper.checkFile(str(p))==str(p)
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(argparse.ArgumentTypeError):
        helper.checkFile(str(p))

def test_closeFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    ps = [d.join('test'+str(ii)+'.txt') for ii in range(10)]
    handles = [helper.openNormalOrGz(str(p),'w') for p in ps]
    helper.closeFiles(handles)
    for ii in handles:
        with pytest.raises(ValueError):
            ii.write()
    handles = dict(zip(range(10),[helper.openNormalOrGz(str(p),'w') for p in ps]))
    helper.closeFiles(handles)
    for _,ii in handles.items():
        with pytest.raises(ValueError):
            ii.write()

def test_checkDir(tmpdir):
    d = tmpdir.mkdir('dir')
    assert helper.checkDir(str(d))==str(d)
    p = d.join('test.txt')
    with pytest.raises(argparse.ArgumentTypeError):
        helper.checkDir(str(p),False)
    p.write("test")
    with pytest.raises(argparse.ArgumentTypeError):
        helper.checkDir(str(p))
    #make unwriteable
    os.chmod(str(d),os.stat(str(d)).st_mode & ~stat.S_IWRITE)
    with pytest.raises(argparse.ArgumentTypeError):
        helper.checkDir(str(d))
    d2=os.path.join(str(tmpdir),'test_make_dir')
    assert not os.path.exists(str(d2))
    assert helper.checkDir(str(d2))==str(d2)
    assert os.path.exists(str(d2))

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


