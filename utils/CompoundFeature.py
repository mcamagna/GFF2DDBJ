from utils.Feature import Feature
"""A special feature that consists of several sub-features. Primarily, this is used for CDS sequences, which can span over multiple introns"""
class CompoundFeature(Feature):

    def __init__(self, members):
        self.members = list(members)
        self.members.sort(key=lambda x: x.start)
        
        start = min([m.start for m in self.members])
        end = max([m.end for m in self.members])
        
        Feature.__init__(self, self.members[0].seqid, self.members[0].source, self.members[0].gfftype, start, end, self.members[0].score, self.members[0].strand, self.members[0].phase, self.members[0].attributes)
        
        self.parent = self.members[0].parent
        children = set()
        for m in self.members:
            children.update(m.children)
        self.children = list(children)
    
    
    def buildLocationString(self):
        """Uses the start/end/strand values to build a location string"""
        if len(self.members) == 1:
            return self.members[0].buildLocationString()
        else:
            s=""
            for m in self.members:
                ms = m.buildLocationString() #need to use the location string, not m.start..m.end since the feature could be truncated
                #only keep the content inside the bracket
                if '(' in ms:
                    ms = ms.split("(", maxsplit=1)[1]
                    ms = ms.split(")", maxsplit=1)[0]
                #s+= str(m.start) + '..'+str(m.end) + ","
                s+= ms+","
            s = s[0:-1] #remove last comma
            s = "join("+s+")"
            if self.strand == '-':
                s= "complement("+s+")"
            return s