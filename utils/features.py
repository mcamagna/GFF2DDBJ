import copy
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
        self.phase = str(phase)
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
    
    
    def _calculatePhase(self):
        if not self.phase.isdigit():
            self.phase = str(0)
            
        
    def split(self, splitstart, splitend):
        """This will split the feature into two pieces, and return two truncated features."""
        left = TruncatedRightFeature.cloneFeature(self)
        right = TruncatedLeftFeature.cloneFeature(self)
        left.end = splitstart-1
        right.start = splitend+1
        left._calculatePhase()
        left.attributes["codon_start"] = str(int(left.phase)+1)
        right._calculatePhase()
        right.attributes["codon_start"] = str(int(right.phase)+1)
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
    
    def removeDownstreamChild(self, obj):
        """Removes a feature/object from the children list of either this object, 
        or the downstream child that is the parent of this object"""
        try:
            self.children.remove(obj)
            return True
        except:
            for child in self.children:
                success = child.removeDownstreamChild(obj)
                if success:
                    return True
        return False
        
        
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
        f.phase = str(self.phase)
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
        
        attributes = copy.deepcopy(self.members[0].attributes) 
        Feature.__init__(self, self.members[0].seqid, self.members[0].source, self.members[0].gfftype, start, end, self.members[0].score, self.members[0].strand, self.members[0].phase, attributes)
        
        self.parent = self.members[0].parent
        children = set()

        for m in self.members:
            children.update(m.children)
        self.children = list(children)    
        self._calculatePhase()
    
    
    def _calculatePhase(self):
            
        if len(self.members)==1:
            self.phase = self.members[0].phase
            return
        
        if self.strand == '+':
            if not isinstance(self.members[0], TruncatedLeftFeature) and not isinstance(self.members[0], TruncatedBothSidesFeature):
                self.phase = 0
            else:
                #do the complicated calculation from the end
                spliced_length = 0
                for m in self.members:
                    spliced_length+= (m.end-m.start)+1
                
                while spliced_length>2:
                    spliced_length -= 3
                self.phase = str(spliced_length)
                     
        elif self.strand == '-':
            if not isinstance(self.members[-1], TruncatedRightFeature) and not isinstance(self.members[-1], TruncatedBothSidesFeature):
                self.phase = self.members[-1].phase
            else:
                #do complicated calculation from the start (the stop codon)
                spliced_length = 0
                for m in self.members:
                    spliced_length+= (m.end-m.start)+1
                
                while spliced_length>2:
                    spliced_length -= 3
                self.phase = str(spliced_length)

        else:
            self.phase = str(0)    
            
    def clone(self):
        f = CompoundFeature(self.members)
        f.attributes = copy.deepcopy(self.attributes)
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
        
        left._calculatePhase()
        left.attributes["codon_start"] = str(int(left.phase)+1)
        right._calculatePhase()
        right.attributes["codon_start"] = str(int(right.phase)+1)
        return [left, right]
    
    def containsTruncatedMember(self):
        for m in self.members:
            if isinstance(m, TruncatedFeature):
                return True
        return False
        
    
class TruncatedFeature(Feature):
    """An empty class, serving as interface for all features that are truncated."""
    def __init__(self, seqid="", source="", gfftype="", start=None, end=None, score=None, strand="", phase="", attribute_dict=None):
        Feature.__init__(self, seqid=seqid, source=source, gfftype=gfftype, start=start, end=end, score=score, strand=strand, phase=phase, attribute_dict=attribute_dict)
    
    
    
"""A special feature that, where the actual start position is smaller than the provided start position"""
class TruncatedLeftFeature(TruncatedFeature):
    
    @staticmethod
    def cloneFeature(basefeature):
        """Make a new TruncatedLeftFeature object from a given basefeature"""
        newfeature = TruncatedLeftFeature()
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
        newfeature._calculatePhase()
        return newfeature
    
    def __init__(self):
        TruncatedFeature.__init__(self, seqid="", source="", gfftype="", start=None, end=None, score=None, strand="", phase="", attribute_dict=None)
    
    
    def _calculatePhase(self):
        
        if self.strand == '-':
            #the start codon is instact
            if not self.phase.isdigit():
                self.phase = str(0)
            
        else:
        #the start has been trimmed, so we need to try and obtain the phase from the stop codon at the end
        #note: if this is part of a compound feature, then the compound feature will have to calculate the 
            length = (self.end - self.start)+1 #+1 because: [1,2,3] -> 3-1 = 2; but the length is actually 3
            while length>2:
                length-=3
            self.phase = str(length)
            
        
    def clone(self):
        f = TruncatedLeftFeature()
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
        right = TruncatedLeftFeature.cloneFeature(self)
        left.end = splitstart-1
        right.start = splitend+1
        #recalculate phase since the positions were changed
        left._calculatePhase()
        left.attributes["codon_start"] = str(int(left.phase)+1)
        right._calculatePhase()
        right.attributes["codon_start"] = str(int(right.phase)+1)
        return [left, right]
    
    




class TruncatedRightFeature(TruncatedFeature):
 
    @staticmethod
    def cloneFeature(basefeature):
        """Make a new TruncatedRightFeature object from a given basefeature"""
        newfeature = TruncatedRightFeature()
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
        newfeature._calculatePhase()
        return newfeature
    
    
    def __init__(self):
        TruncatedFeature.__init__(self, seqid="", source="", gfftype="", start=None, end=None, score=None, strand="", phase="", attribute_dict=None)
        self._calculatePhase()
    
    
    def _calculatePhase(self):
            
        if self.strand == '-':
            #the start has been trimmed, so we need to try and obtain the phase from the stop codon at the end
            #note: if this is part of a compound feature, then the compound feature will have to calculate the 
            length = (self.end - self.start)+1 #+1 because: [1,2,3] -> 3-1 = 2; but the length is actually 3
            while length>2:
                length-=3
            
            self.phase = str(length)
        else:
            #the start codon is instact
            if not self.phase.isdigit():
                self.phase = str(0)
    
    
    def clone(self):
        f = TruncatedRightFeature()
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
        left = TruncatedRightFeature.cloneFeature(self)
        right = TruncatedBothSidesFeature.cloneFeature(self)
        left.end = splitstart-1
        right.start = splitend+1
        #recalculate phase since the positions were changed
        left._calculatePhase()
        left.attributes["codon_start"] = str(int(left.phase)+1)
        right._calculatePhase()
        right.attributes["codon_start"] = str(int(right.phase)+1)
        return [left, right]
    
    
class TruncatedBothSidesFeature(TruncatedFeature):
 
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
        newfeature._calculatePhase()
        return newfeature
    
    
    def __init__(self):
        TruncatedFeature.__init__(self, seqid="", source="", gfftype="", start=None, end=None, score=None, strand="", phase="", attribute_dict=None)
        self._calculatePhase()
        
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
    
    def _calculatePhase(self):
        #This is not correct, though there is no way to calculate the phase
        #If start and stopcodons are missing, then we simply can't know the phase if it wasn't provided in the GFF
        if not self.phase.isdigit():
            self.phase = str(0)
    
    def split(self, splitstart, splitend):
        """This will split the feature into two pieces, and return two truncated features."""
        left = TruncatedBothSidesFeature.cloneFeature(self)
        right = TruncatedBothSidesFeature.cloneFeature(self)
        left.end = splitstart-1
        right.start = splitend+1
        #recalculate phase since the positions were changed
        left._calculatePhase()
        left.attributes["codon_start"] = str(int(left.phase)+1)
        right._calculatePhase()
        right.attributes["codon_start"] = str(int(right.phase)+1)
        
        return [left, right]