import pytest
from dnapy import countbases
import argparse
import os, stat

def test_checkFile(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    with pytest.raises(argparse.ArgumentTypeError):
        countbases.check_file(str(d))
    #doesn't exist yet
    with pytest.raises(argparse.ArgumentTypeError):
        countbases.check_file(str(p))
    p.write("test")
    assert countbases.check_file(str(p))==str(p)
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(argparse.ArgumentTypeError):
        countbases.check_file(str(p))


