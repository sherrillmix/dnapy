#based on https://biopython.org/wiki/ABI_traces
#lapply(ab1,function(xx){tmp<-tempfile();system(sprintf('abitotrace %s>%s',xx,tmp));read.csv(tmp)})
#cols=c('G'='black','T'='red','A'='green','C'='blue');pdf('traces.pdf',width=20,height=8);lapply(traces,function(xx){plot(1,1,type='n',xlim=c(1,nrow(xx)),ylim=range(xx),xlab='Time',ylab='Intensity',las=1);for(ii in 1:ncol(xx))lines(xx[,ii],col=cols[colnames(xx)[ii]])});dev.off()
from Bio import SeqIO
from dnapy import helper
import argparse

def main(argv=None):
    parser = argparse.ArgumentParser(description="A program to take an ABI .ab1 file from Sanger sequencing and output the traces to standard out as four columns.")
    parser.add_argument('abiFile', help='an ABI .ab1 file containing the sequence read',type=helper.checkFile)

    args=parser.parse_args(argv)
    record = SeqIO.read(args.abiFile, 'abi')
    channels = {'DATA9': 'G', 'DATA10':'A', 'DATA11':'T', 'DATA12':'C'}
    
    traces={channels[chan]:[str(ii) for ii in record.annotations['abif_raw'][chan]] for chan in channels.keys()}
    lengths=[len(traces[xx]) for xx in traces]
    if(any([l!=lengths[1] for l in  lengths])): raise IndexError('Length of traces inconsistent')
    out=[','.join(gatc) for gatc in zip(traces['G'],traces['A'],traces['T'],traces['C'])]

    print("G,A,T,C")
    for ii in out: print(ii)
    

if __name__ == '__main__':
    main()
