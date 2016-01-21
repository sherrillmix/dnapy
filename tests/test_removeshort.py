import pytest
from dnapy import removeshort
import os
import stat
import pysam
import subprocess

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    with pytest.raises(IOError):
        next(removeshort.removeShort(str(d)))
    #doesn't exist yet
    with pytest.raises(IOError):
        next(removeshort.removeShort(str(p)))
    #incorrectly formatted file
    #doesn't exist yet
    p.write("test")
    #with pytest.raises(ValueError):
    #removeshort.removeShort(str(p))
    #doesn't exist yet
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(IOError):
        next(removeshort.removeShort(str(p)))

