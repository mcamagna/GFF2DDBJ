"""The Feature class stores all required information for features, including their relationships."""  
class Feature:
  
    def __init__(self, seqid="", source="", gfftype="", start=None, end=None, score=None, strand="", phase="", attribute_dict=None):
        self.seqid = seqid
        self.source = source
        self.gfftype = gfftype
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = attribute_dict
        if self.attributes is None:
            self.attributes = dict()
        
        self.parent = None #reference to the parent feature object
        self.children = [] #list of feature objects belonging to this feature
    
    def buildLocationString(self):
        """Uses the start/end/strand values to build a location string"""
        if self.start == self.end:
            return str(self.start)
        
        s = str(self.start) + '..'+str(self.end)
        if self.strand == '-':
            s= "complement("+s+")"
        return s
    
    def split(self, splitstart, splitend):
        """This will split the feature into two pieces, and return two truncated features."""
        left = TruncatedEndFeature.cloneFeature(self)
        right = TruncatedStartFeature.cloneFeature(self)
        left.end = splitstart-1
        right.start = splitend+1
        return [left, right]
    
    def addAttribute(self, name, value):
        self.attributes[name] = str(value)
    
    def getAttribute(self, name):
        return self.attributes.get(name)
    
    def hasAttribute(self, name, case_sensitive=False):
        """Checks whether the feature has a given attribute."""
        if not case_sensitive:
            name = name.upper()
            return name.upper() in [k.upper() for k in self.attributes.keys()]
        else:
            return k in self.attributes.keys()
        
            
    def getParent(self):
        return self.parent
    
    def getUltimateParent(self):
        if self.parent is not None:
            return self.parent.getUltimateParent()
        return self
    
    def getAllDownstreamChildren(self):
        allchildren = set()
        for child in self.children:
            allchildren.add(child)
            allchildren.update(child.getAllDownstreamChildren())
        return allchildren
    
    def getAllDownstreamOfType(self, gfftypes):
        if not isinstance(gfftypes, list) and not isinstance(gfftypes, set):
            gfftypes = [gfftypes]
            
        to_return = set()
        for child in self.children:
            if child.gfftype in gfftypes:
                to_return.add(child)
                to_return.update(child.getAllDownstreamOfType(gfftypes))
        return to_return
        
    def getAllDownstreamCDS(self):
        return self.getAllDownstreamOfType(["CDS", "cds"])
    
    def getAllDownstreamMRNA(self):
        return self.getAllDownstreamOfType(["mRNA", "MRNA", "mrna"])
    
    def sortChildrenByPosition(self):
        self.children.sort(key=lambda x: (x.start, x.end))
    
    def getHash(self):
        h = "" + self.seqid + '_' + self.source + '_' + self.gfftype 
        h+= '_' + self.buildLocationString()
        h+= '_' +self.strand
        attr_items = list(self.attributes.items())
        attr_items.sort(key=lambda x: x[0])
        for qualifier, value in attr_items:
            h+="_"+str(qualifier)+"_"+str(value)
        return h
    
    def removeDuplicateChildren(self):
        children_to_remove = []
        child_hashes = set()
        for child in self.children:
            h = child.getHash()
            
            if h in child_hashes:
                children_to_remove.append(child)
            else:
                child_hashes.add(h)
        for ctr in children_to_remove:
            self.children.remove(ctr)
        
    
    
    
    def clone(self):
        f = Feature()
        import copy
        f.seqid = self.seqid
        f.source = self.source
        f.gfftype = self.gfftype
        f.start = self.start
        f.end = self.end
        f.score = self.score
        f.strand = self.strand
        f.phase = self.phase
        f.attributes = copy.deepcopy(self.attributes)    
        f.parent = self.parent
        f.children = self.children
        return f
    
    
    
    
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
    
    def clone(self):
        f = CompoundFeature(self.members)
        return f
    
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
        
        
    def split(self, splitstart, splitend):
        """This will split the feature into two pieces, and return two truncated features."""
        if len(self.members)==1:
            splitmember = self.members[0].split(splitstart, splitend)
            return [splitmember[0], splitmember[1]]
        
        left_members = []
        right_members = []
        #note: members are already sorted by position
        for member in self.members:
            if member.end<splitstart:
                left_members.append(member)
            elif member.start>splitend:
                right_members.append(member)
            else:
                splitmember = member.split(splitstart, splitend)
                left_members.append(splitmember[0])
                right_members.append(splitmember[1])
        
        left = self.clone()
        left.members = left_members
        right = self.clone()
        right.members = right_members
        left.end = splitstart-1
        right.start = splitend+1
        return [left, right]
    
    
    
    
"""A special feature that, where the actual start position is smaller than the provided start position"""
class TruncatedStartFeature(Feature):
    
    @staticmethod
    def cloneFeature(basefeature):
        """Make a new TruncatedStartFeature object from a given basefeature"""
        newfeature = TruncatedStartFeature()
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
    
    
    def clone(self):
        f = TruncatedStartFeature()
        import copy
        f.seqid = self.seqid
        f.source = self.source
        f.gfftype = self.gfftype
        f.start = self.start
        f.end = self.end
        f.score = self.score
        f.strand = self.strand
        f.phase = self.phase
        f.attributes = copy.deepcopy(self.attributes)    
        f.parent = self.parent
        f.children = self.children
        return f
    
    
    def buildLocationString(self):
        s = "<"+str(self.start) + '..'+str(self.end)
        if self.strand == '-':
                s = "complement("+s+")"
        return s
    
    
    def split(self, splitstart, splitend):
        """This will split the feature into two pieces, and return two truncated features."""
        left = TruncatedBothSidesFeature.cloneFeature(self)
        right = TruncatedStartFeature.cloneFeature(self)
        left.end = splitstart-1
        right.start = splitend+1
        return [left, right]
    
    




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
    
    def clone(self):
        f = TruncatedEndFeature()
        import copy
        f.seqid = self.seqid
        f.source = self.source
        f.gfftype = self.gfftype
        f.start = self.start
        f.end = self.end
        f.score = self.score
        f.strand = self.strand
        f.phase = self.phase
        f.attributes = copy.deepcopy(self.attributes)    
        f.parent = self.parent
        f.children = self.children
        return f
    
    
    def buildLocationString(self):
        s = str(self.start) + '..'+">"+str(self.end)
        if self.strand == '-':
            s = "complement("+s+")"
        return s
    
    
    def split(self, splitstart, splitend):
        """This will split the feature into two pieces, and return two truncated features."""
        left = TruncatedEndFeature.cloneFeature(self)
        right = TruncatedBothSidesFeature.cloneFeature(self)
        left.end = splitstart-1
        right.start = splitend+1
        return [left, right]
    
    
class TruncatedBothSidesFeature(Feature):
 
    @staticmethod
    def cloneFeature(basefeature):
        """Make a new TruncatedBothSidesFeature object from a given basefeature"""
        newfeature = TruncatedBothSidesFeature()
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
    
    def clone(self):
        f = TruncatedBothSidesFeature()
        import copy
        f.seqid = self.seqid
        f.source = self.source
        f.gfftype = self.gfftype
        f.start = self.start
        f.end = self.end
        f.score = self.score
        f.strand = self.strand
        f.phase = self.phase
        f.attributes = copy.deepcopy(self.attributes)    
        f.parent = self.parent
        f.children = self.children
        return f
    
    
    def buildLocationString(self):
        s = '<'+str(self.start) + '..'+">"+str(self.end)
        if self.strand == '-':
            s = "complement("+s+")"
        return s
    
    
    def split(self, splitstart, splitend):
        """This will split the feature into two pieces, and return two truncated features."""
        left = TruncatedBothSidesFeature.cloneFeature(self)
        right = TruncatedBothSidesFeature.cloneFeature(self)
        left.end = splitstart-1
        right.start = splitend+1
        return [left, right]