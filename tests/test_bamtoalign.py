import pytest
from dnapy import helper
from dnapy import bamtoalign
import pysam


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
    with pytest.raises(ImportError):
        bamtoalign.getRefFromFasta([['Test1','AAAACCCCTTTTGGGG'],['Test2','ACACA']],'Test')
    with pytest.raises(ImportError):
        bamtoalign.getRefFromFasta([['Test','AAAACCCCTTTTGGGG'],['Test','ACACA']],'Test')
    seq=bamtoalign.getRefFromFasta([['Test','AAAACCCCTTTTGGGG'],['Test2','ACACA']],'Test')
    assert seq[0]=='Test'
    assert seq[1]=='AAAACCCCTTTTGGGG'
    seq=bamtoalign.getRefFromFasta([['Test','AAAACCCCTTTTGGGG'],['Test2','ACACA']],'Test2')
    assert seq[0]=='Test2'
    assert seq[1]=='ACACA'

#something about capsys messes up pysam.index so create ahead of time for testing main
@pytest.fixture(scope='session')
def bamFile(tmpdir_factory):
    header = { 'HD': {'VN': '1.0'}, 'SQ': [{'LN': 1000, 'SN': 'ref'}] }
    p=tmpdir_factory.mktemp('test').join('test.bam')
    outFile=pysam.AlignmentFile(str(p),"wb",header=header)
    a = pysam.AlignedSegment()
    a.query_name = "read3"
    a.query_sequence="GGGGAAAAAT"
    a.reference_start = 28
    a.reference_id = 0
    a.mapping_quality = 20
    a.cigar = ((0,10), )
    a.flag=16
    outFile.write(a)
    a.query_name = "read2"
    a.reference_start = 32
    a.query_sequence="AAAAATTTTT"
    a.cigar = ((1,1),(2,2),(0,8),(1,1))
    a.flag=0
    outFile.write(a)
    a.mapping_quality = 19
    a.query_name = "read1"
    a.query_sequence="TTCCAAAAACCCCCGGCC"
    a.cigar = ((5,5),(4,2),(1,2),(0,9),(8,1), (2,1),(3,1),(7,1),(1,1),(0,2))
    outFile.write(a)
    outFile.close()
    pysam.index(str(p))
    return(p)


def test_trimCigar():
    assert bamtoalign.trimCigar([[0,10]],0,3)==[[[0,10]],3]
    assert bamtoalign.trimCigar([[0,1]],0,5)==[[[0,1]],5]
    assert bamtoalign.trimCigar([[1,2],[0,10],[1,2]],0,3)==[[[4,2],[0,10],[4,2]],3]
    assert bamtoalign.trimCigar([[2,2],[0,10],[2,2]],0,3)==[[[6,2],[0,10],[6,2]],5]
    assert bamtoalign.trimCigar([[0,10],[2,2]],0,3)==[[[0,10],[6,2]],3]
    assert bamtoalign.trimCigar([[0,2],[2,2],[0,10],[2,2]],2,3)==[[[4,2],[6,2],[0,10],[6,2]],7]
    assert bamtoalign.trimCigar([[0,2],[2,2],[0,10],[2,2]],1,3)==[[[0,2],[2,2],[0,10],[6,2]],3]

def test_getAlignsInFile(tmpdir,bamFile):
    for xx,yy in zip(bamtoalign.getAlignsInFile(str(bamFile)),[
        {'name':'read3','start':28,'seq':'GGGGAAAAAT','strand':'-','insertions':[]},
        {'name': 'read2', 'start': 32, 'insertions': [[32,'A'],[42,'T']], 'strand': '+', 'seq': '--AAAATTTT'},
        {'name': 'read1', 'start': 32, 'insertions': [[32,'CC'],[45,'G']], 'strand': '+', 'seq': 'AAAAACCCCC--GCC'}
    ]):
        assert xx==yy
    for xx,yy in zip(bamtoalign.getAlignsInFile(str(bamFile),endSpan=1),[
        {'name':'read3','start':28,'seq':'GGGGAAAAAT','strand':'-','insertions':[]},
        {'name': 'read2', 'start': 32, 'insertions': [], 'strand': '+', 'seq': 'AAAATTTT'},
        {'name': 'read1', 'start': 32, 'insertions': [[45,'G']], 'strand': '+', 'seq': 'AAAAACCCCC--GCC'}
    ]):
        assert xx==yy
    for xx,yy in zip(bamtoalign.getAlignsInFile(str(bamFile),endSpan=2),[
        {'name':'read3','start':28,'seq':'GGGGAAAAAT','strand':'-','insertions':[]},
        {'name': 'read2', 'start': 32, 'insertions': [], 'strand': '+', 'seq': 'AAAATTTT'},
        {'name': 'read1', 'start': 32, 'insertions': [], 'strand': '+', 'seq': 'AAAAACCCC'}
    ]):
        assert xx==yy

    for xx,yy in zip(bamtoalign.getAlignsInFile(str(bamFile),'ref:28'),[{'name':'read3','start':28,'seq':'GGGGAAAAAT','strand':'-','insertions':[]}]):
        assert xx==yy
    for xx,yy in zip(bamtoalign.getAlignsInFile(str(bamFile),minQuality=20),[
        {'name':'read3','start':28,'seq':'GGGGAAAAAT','strand':'-','insertions':[]},
        {'name': 'read2', 'start': 32, 'insertions': [[32,'A'],[42,'T']], 'strand': '+', 'seq': '--AAAATTTT'}
    ]):
        assert xx==yy
    assert len([xx for xx in bamtoalign.getAlignsInFile(str(bamFile),'ref:900')])==0
    assert len([xx for xx in bamtoalign.getAlignsInFile(str(bamFile),minQuality=21)])==0
    with pytest.raises(ValueError):
        [xx for xx in bamtoalign.getAlignsInFile(str(bamFile),'notRealRef:1')]


def test_main(capsys,tmpdir,bamFile):
    with pytest.raises(SystemExit):
        bamtoalign.main()
    out, err=capsys.readouterr()
    assert 'usage' in err
    with pytest.raises(SystemExit):
        bamtoalign.main(['-h'])
    out, err=capsys.readouterr()
    assert 'usage' in out
    d = tmpdir.mkdir('dir')
    p = d.join('ref.fa')
    with helper.openNormalOrGz(str(p),'w') as f:
        f.write(">ref\nAAAAAAAAAATTTTTTTTTTCCCCCCCCCCGGGGGGGGGGAAAAAAAAAA")
    bamtoalign.main(['-v','-s',str(p),str(bamFile)])
    out, err=capsys.readouterr()
    assert 'Arguments' in err
    for ii,jj in zip(out.split('\n'),[
        '>ref',  'AAAAAAAAAATTTTTTTTTTCCCCCCCCCCGG--GGGGGGGGAA-AAA-AAAAA',
        '>read3','----------------------------GGGG--AAAAAT--------------',
        '>read2','--------------------------------A---AAAATTTTT---------',
        '>read1','--------------------------------CCAAAAACCCCC---GGCC---'
    ]):
        assert ii==jj


