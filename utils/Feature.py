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
        self.children.sort(key=lambda x: x.start)
    
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