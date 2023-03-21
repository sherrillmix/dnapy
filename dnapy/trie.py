class Node:
    def __init__(self,parent):
        self.children={}
        self.end=False

class Trie:
    def __init__(self):
        self.root=Node(None)

    def insert(self,seq):
        node=self.root
        for char in seq:
            if char in node.children:
                node=node.children[char]
            else:
                node.children[char]=Node(node)
                node=node.children[char]
        node.end=True
    
    def getSeqs(self,node=None):
        if node is None:
            node=self.root
        if len(node.children)==0:
            return ['']
        else:
            return [child+seq for child in node.children for seq in self.getSeqs(node.children[child])]

    def checkSeq(self,seq,node=None):
        if node is None: node=self.root
        for char in seq:
            if char in node.children:
                node=node.children[char]
            else:
                return False
        return node.end

    def __contains__(self, key):
        return self.checkSeq(key)

    def checkError(self,seq,node=None,maxErrors=1):
        if node is None: node=self.root
        if seq=='':
            if node.end:
                return [('',0)]
            else:
                return []
        if maxErrors<1:
            if self.checkSeq(seq,node):
                return [(seq,0)]
            else:
                return []
        else:
            #for each possible child node check downstream and if possible return tuple of (childChar + possibleString, numberErrorsInPossibleString) if childChar matches seq[0] else (childChar + possibleString, numberErrorsInPossibleString+1)
            out=[(child+ii[0],ii[1] if child == seq[0] else ii[1]+1) for child in node.children for ii in self.checkError(seq[1:],node.children[child],maxErrors if child == seq[0] else  maxErrors-1) if ii]
            if len(out)>0:
                return out
            else:
                return []


