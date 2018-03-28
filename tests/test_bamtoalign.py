import pytest
from dnapy import helper
from dnapy import bamtoalign


def test_padRead():
    assert bamtoalign.padRead("ACACA",0,5)=='ACACA'
    with pytest.raises(IndexError):
        bamtoalign.padRead("ACACA",1,5)
    assert bamtoalign.padRead("ACACA",1,7)=='-ACACA-'
    assert bamtoalign.padRead("AC-CA",2,7)=='--AC-CA'
    assert bamtoalign.padRead("AC-CA",2,7)=='--AC-CA'
    assert bamtoalign.padRead("ACACA",0,5,{},{1:2})=='A--CACA'
    assert bamtoalign.padRead("ACACA",0,6,refInserts={1:4})=='A----CACA-'
    assert bamtoalign.padRead("ACACA",0,5,[(1,'ZZ')],{1:4})=='AZZ--CACA'
    with pytest.raises(IndexError):
        bamtoalign.padRead("ACACA",0,5,[(2,'ZZ')],{1:4})
    assert bamtoalign.padRead("ACACA",0,5,[(2,'ZZ')],{2:2})=='ACZZACA'
    assert bamtoalign.padRead("ACACA",3,8,[(2,'ZZ')],{2:2})=='--ZZ-ACACA'
    with pytest.raises(IndexError):
        bamtoalign.padRead("ACACA",0,5,[(2,'ZZ')],{2:1})

def test_getRefFromFasta():
    with pytest.raises(ImportError):
        bamtoalign.getRefFromFasta([])
    seq=bamtoalign.getRefFromFasta([['Test','AAAACCCCTTTTGGGG']])
    assert seq[0]=='Test'
    assert seq[1]=='AAAACCCCTTTTGGGG'
    with pytest.raises(ImportError):
        bamtoalign.getRefFromFasta([['Test','AAAACCCCTTTTGGGG'],['Test2','ACACA']])
    seq=bamtoalign.getRefFromFasta([['Test','AAAACCCCTTTTGGGG'],['Test2','ACACA']],'Test')
    assert seq[0]=='Test'
    assert seq[1]=='AAAACCCCTTTTGGGG'
    seq=bamtoalign.getRefFromFasta([['Test','AAAACCCCTTTTGGGG'],['Test2','ACACA']],'Test2')
    assert seq[0]=='Test2'
    assert seq[1]=='ACACA'



