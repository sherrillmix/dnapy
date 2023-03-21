import pytest
from dnapy import trie

def test_trie():
    #with pytest.raises(argparse.ArgumentTypeError):
        #helper.checkPositiveInt(10,9)

    t=trie.Trie()
    assert not t.checkSeq('ABD')
    inputSeqs=sorted(['AB1','ACC'])
    for ii in inputSeqs: t.insert(ii)
    seqs=sorted(t.getSeqs())
    assert seqs==inputSeqs
    assert t.checkSeq('AB1')
    assert t.checkSeq('ACC')
    assert not t.checkSeq('ABD')
    assert not t.checkSeq('AB11')
    assert not t.checkSeq('AB1C')
    assert not t.checkSeq('AB')
    assert not t.checkSeq('')
    assert t.checkError('')==[]
    assert t.checkError('ACC1')==[]
    assert t.checkError('ABD',maxErrors=0)==[]
    assert sorted(t.checkError('ABD',maxErrors=1))==sorted([('AB1',1)])
    assert sorted(t.checkError('ABD',maxErrors=2))==sorted([('AB1',1),('ACC',2)])
    assert sorted(t.checkError('AB1',maxErrors=2))==sorted([('AB1',0),('ACC',2)])
    t.insert('ACC1')
    assert t.checkError('ACC1')==[('ACC1',0)]
    assert t.checkError('ACCZ',maxErrors=2)==[('ACC1',1)]
    assert sorted(t.checkError('AB1',maxErrors=2))==sorted([('AB1',0),('ACC',2)])


