import pytest
from dnapy import splitbarcodes
from dnapy import helper
import os
import stat

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p = d.join('test.txt')
    p1 = d.join('test1.txt')
    p1.write("@test\nAAAA\n+\n1111")
    #not a file
    with pytest.raises(IOError):
        next(splitbarcodes.barcodeFastqIter([str(d)],[str(p1)],set()))
    with pytest.raises(IOError):
        next(splitbarcodes.barcodeFastqIter([str(p1)],[str(d)],set()))
    #doesn't exist yet
    with pytest.raises(IOError):
        splitbarcodes.barcodeFastqIter([str(p)],[str(p1)],set())
    with pytest.raises(IOError):
        splitbarcodes.barcodeFastqIter([str(p1)],[str(p)],set())
    #incorrectly formatted file
    p.write("@test\nAAAA\n+\n1") #qual and seq different
    with pytest.raises(ValueError):
        next(splitbarcodes.barcodeFastqIter([str(p)],[str(p1)],set()))
    with pytest.raises(ValueError):
        next(splitbarcodes.barcodeFastqIter([str(p1)],[str(p)],set()))
    #make unreadable
    os.chmod(str(p),os.stat(str(p)).st_mode & ~stat.S_IREAD)
    with pytest.raises(IOError):
        splitbarcodes.barcodeFastqIter([str(p)],[str(p1)],set())
    with pytest.raises(IOError):
        splitbarcodes.barcodeFastqIter([str(p1)],[str(p)],set())
