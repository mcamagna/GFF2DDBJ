from utils.Feature import Feature
"""A special feature that, where the actual end position is higher than the provided end position"""
class TruncatedEndFeature(Feature):
 
    @staticmethod
    def cloneFeature(basefeature):
        """Make a new TruncatedEndFeature object from a given basefeature"""
        newfeature = TruncatedEndFeature()
        newfeature.seqid = basefeature.seqid
        newfeature.source = basefeature.source
        newfeature.gfftype = basefeature.gfftype
        newfeature.start = basefeature.start
        newfeature.end = basefeature.end
        newfeature.score = basefeature.score
        newfeature.strand = basefeature.strand
        newfeature.phase = basefeature.phase
        newfeature.attributes = basefeature.attributes.copy() 
        newfeature.parent = basefeature.parent
        newfeature.children = basefeature.children #list of feature objects belonging to this feature
        return newfeature
    
    
    def __init__(self):
        Feature.__init__(self, seqid="", source="", gfftype="", start=None, end=None, score=None, strand="", phase="", attribute_dict=None)
    
    def buildLocationString(self):
        s = str(self.start) + '..'+">"+str(self.end)
        if self.strand == '-':
            s = "complement("+s+")"
        return s