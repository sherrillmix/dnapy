#abi file from https://github.com/bow/abifpy/blob/master/tests/310.ab1
import pytest
from dnapy import helper
from dnapy import abitotrace
import os

def test_badFiles(tmpdir):
    d = tmpdir.mkdir('dir')
    p1 = d.join('test.txt')
    #not a file
    with pytest.raises(SystemExit):
        abitotrace.main([str(d)])
    #doesn't exist yet
    with pytest.raises(SystemExit):
        abitotrace.main([str(p1)])
    #incorrectly formatted file
    p1.write("1234567890") #qual and seq different
    with pytest.raises(IOError):
        abitotrace.main([str(p1)])

def test_realAbi(capsys):
    ab1=os.path.join(os.path.dirname(os.path.abspath(__file__)),'310.ab1')
    abitotrace.main([ab1])
    out, err=capsys.readouterr()
    assert err==''
    outSplit=out.split('\n')
    assert outSplit[0]=='G,A,T,C'
    assert outSplit[1]=='0,0,0,115'
    assert outSplit[3]=='0,0,0,73'
    assert outSplit[9826]=='5,10,13,6'
    assert outSplit[1000]=='126,40,0,49'

