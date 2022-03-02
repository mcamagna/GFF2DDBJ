'''
Braker2 generated GFF files lack the gene feature. We'll try to parse the file as generic as possible to allow conversion of other GFF files as well
@author: Maurizio Camagna
'''

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
        self.attributes = None
        if self.attributes is None:
            self.attributes = dict()
        
        self.parent = None #reference to the parent feature object
        self.children = [] #list of feature objects belonging to this feature
    
    def addAttribute(self, name, value):
        self.attributes[name] = str(value)
    
    def getAttribute(self, name):
        return self.attributes.get(name)
    
    
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
    
    
class GFFParser:
    
    def __init__(self, gff_path):
        self.gff_path = gff_path
        self.features = {} #contains all features
        self.parentFeatures = [] #contains only features that have no parent themselves
    
        self._parseGFF()
        self._connectFeatures()
    
    
    def _parseGFF(self):
        """Iterates through a GFF file and extracts all features."""
        file_handle = None
        if self.gff_path.lower().endswith(".gz"):
            import gzip
            file_handle = gzip.open(self.gff_path, 'rt')
        else:
            file_handle = open(self.gff_path, 'rt')
        
        for line in file_handle:
            if line.startswith("#"):
                continue
            line = line.replace("\n", "")
            #print(line)
            spl = line.split("\t")
            seqid=spl[0]
            source=spl[1]
            gfftype=spl[2]
            start=int(spl[3])
            end=int(spl[4])
            score=spl[5]
            strand=spl[6]
            phase=spl[7]
            attributes = spl[8].split(";")
            
            feature = Feature(seqid=seqid, source=source, gfftype=gfftype, start=start, end=end, score=score, strand=strand,phase=phase)
            
            for attribute in attributes:
                if "=" in attribute:
                    attsplit = attribute.split("=")
                    name = attsplit[0]
                    value = attsplit[1]
                    feature.addAttribute(name, value)
            
            #print(feature.getAttribute("ID"), feature.getAttribute("Parent"),attributes)

            self.features[feature.getAttribute("ID")] = feature 
            
        file_handle.close()
        
        
    def _connectFeatures(self):    
        """Connects features that are in a parent/child relationship."""
        
        placeholders = [] #a list of placeholder features, which we may need of some parent features are missing
        
        keys = list(self.features.keys()) #save a copy because the size of the dict may change
        
        for f_id in keys:
            #print(f_id)
            f = self.features.get(f_id)
            parent_id = f.getAttribute("Parent")
            if parent_id is None:
                self.parentFeatures.append(f)
            else:
                #check if the parent exists
                parent = self.features.get(parent_id)
                #print(parent_id)
                if parent is None:
                    #some GFF files, including Braker2 generated ones don't actually
                    #provide the parent. We'll have to create the parent in that case
                    #For now we'll make a placeholder feature. Once all features are connect, we can
                    #calculate the start/end/strand etc from all the child nodes
                    parent = Feature(gfftype="PLACEHOLDER")
                    placeholders.append(parent)
                    parent.addAttribute("ID", parent_id)
                    self.features[parent_id] = parent
                    self.parentFeatures.append(parent)
                    
                parent.children.append(f)
                f.parent = parent
        
        for placeholder in placeholders:
            self._fixPlaceholder(placeholder)
                
    
           
    def _fixPlaceholder(self, placeholder_feature):
        """Placeholder feauters are parent features that were missing in the GFF file.
        Using their child nodes, we can infer the required information.
        """
        children_contain_gene_feature = False
        
        child_start_positions = []
        child_end_positions = []
        strand = None
        seqid = None
        source = None
        
        for child in placeholder_feature.children:
            child_start_positions.append(child.start)
            child_end_positions.append(child.end)
            
            child_type = child.gfftype.upper()
            if child_type == "GENE":
                children_contain_gene_feature = True
            
            if seqid is None:
                seqid = child.seqid
            if source is None:
                source = child.source
                
            if strand is None:
                strand = child.strand
            elif (child.strand == "+" or child.strand == "-") and (strand != '+' and strand != "-"):
                #if the child has + or - let's use that over a . entry
                strand = child.strand
            
        placeholder_feature.seqid = seqid
        placeholder_feature.source = source
        placeholder_feature.start = min(child_start_positions)
        placeholder_feature.end = max(child_end_positions)
        placeholder_feature.strand = strand
        placeholder_feature.score = '.'
        placeholder_feature.phase = '.'
        
        
        if not children_contain_gene_feature:
            placeholder_feature.gfftype = "gene"
        else:
            print("Warning: Found gene with an unknown parent. Adding 'source' feature as parent of that gene")
            placeholder_feature.gfftype = "source"
        
           
    