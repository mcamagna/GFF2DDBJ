'''
Braker2 generated GFF files lack the gene feature. We'll try to parse the file as generic as possible to allow conversion of other GFF files as well
@author: Maurizio Camagna
'''
from utils.DDBJWriter import DDBJWriter
from utils.features import Feature, CompoundFeature, TruncatedLeftFeature, TruncatedRightFeature,\
    TruncatedBothSidesFeature
from utils.Parameters import Parameters

    
    
class GFFParser:
    
    def __init__(self, gff_path):
        self.gff_path = gff_path
        self.features = {} #contains all features
        self.parentFeatures = [] #contains only features that have no parent themselves
    
        self._parseGFF()
        self._connectFeatures()
        self._mergeCDSSequences()
        self._detectIncompleteCDS()
    
    
    def _preparseGFF(self, lines=1000):
        """Reads the start of the GFF file to determine what type of GFF3 file is present"""
        file_handle = None
        if self.gff_path.lower().endswith(".gz"):
            import gzip
            file_handle = gzip.open(self.gff_path, 'rt')
        else:
            file_handle = open(self.gff_path, 'rt')
        
        for i, line in enumerate(file_handle):
            if i>lines:
                break
            if "\tstart_codon\t" in line:
                Parameters.gff_contains_startcodons = True
            elif "\tgene\t" in line:
                Parameters.gff_contains_genes = True
            elif "\tmRNA\t" in line:
                Parameters.gff_contains_transcripts = True
            
            #lets check if we already found all information we were looking for and stop if so
            if Parameters.gff_contains_startcodons and Parameters.gff_contains_genes and Parameters.gff_contains_transcripts:
                break
            
        file_handle.close()
    
    
    def _parseGFF(self):
        """Iterates through a GFF file and extracts all features."""
        
        self._preparseGFF() #pre-parse the GFF to see what type of GFF file it is
        
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
            
            
            #Some GFF files assign the same ID to all CDS fragments, spread over multiple lines, others use different ID's
            #Sometimes, the ID is even completely missing for child nodes
            #If the ID is already present, we can assign a new ID and later merge them via their shared parent
            if feature.getAttribute("ID") is None:
                if feature.getAttribute("Parent") is not None:
                    feature.addAttribute("ID", feature.getAttribute("Parent")+"X")
                pass
            while feature.getAttribute("ID") in self.features.keys():
                feature.attributes["ID"] = feature.attributes["ID"]+'X'
            self.features[feature.getAttribute("ID")] = feature
            
        file_handle.close()
        
        
    def _connectFeatures(self):    
        """Connects features that are in a parent/child relationship."""
        
        placeholders = [] #a list of placeholder features, which we may need if some parent features are missing
        
        keys = list(self.features.keys()) #save a copy because the size of the dict may change
        
        for f_id in keys:
            #print(f_id)
            feature = self.features.get(f_id)
            
            parent_id = feature.getAttribute("Parent")
            if parent_id is None:
                self.parentFeatures.append(feature)
            else:
                #check if the parent exists
                parent = self.features.get(parent_id)
                if parent is not None:
                    parent.children.append(feature)
                    feature.parent = parent
                else:
                    #some GFF files, don't actually
                    #provide the parent. We'll have to create the parent in that case
                    #For now we'll make a placeholder feature. Once all features are connected, we can
                    #calculate the start/end/strand etc from all the child nodes
                    if Parameters.gff_contains_transcripts and Parameters.gff_contains_genes == False:
                        #the gff file contains mRNAs, so the missing parent must belong be a gene
                        parent = Feature(gfftype="gene")
                        placeholders.append(parent)
                        parent.addAttribute("ID", parent_id)
                        self.features[parent_id] = parent
                        self.parentFeatures.append(parent)
                        parent.children.append(feature)
                        feature.parent = parent
                        
                        
                    elif Parameters.gff_contains_transcripts == False and Parameters.gff_contains_genes == False:
                        #the gff file contains neither gene nor mRNA. We will need to create both
                        mRNA = Feature(gfftype="mRNA")
                        placeholders.append(mRNA)
                        gene = Feature(gfftype="gene")
                        placeholders.append(gene)
                        
                        mRNA.addAttribute("ID", parent_id)
                        mRNA.children.append(feature)
                        feature.parent = mRNA
                        
                        gene_id = None
                        if "." in parent_id:
                            gene_id = parent_id.rsplit(".", maxsplit=1)[0]
                        else:
                            gene_id = "gene_"+parent_id
                        gene.addAttribute("ID", gene_id)
                        gene.children.append(mRNA)
                        mRNA.parent = gene
                        self.parentFeatures.append(gene)
                        
                        self.features[gene.getAttribute("ID")] = gene
                        self.features[mRNA.getAttribute("ID")] = mRNA
                    
                    
        for placeholder in placeholders:
            self._fixPlaceholder(placeholder)
        
    
           
    def _fixPlaceholder(self, placeholder_feature):
        """Placeholder features are features, that were created because the entries in the GFF
        file were referring to a parent, that was not written into the GFF file.
        Using their child nodes, we can infer the required information.
        """
        
        child_start_positions = []
        child_end_positions = []
        strand = None
        seqid = None
        source = None
        
        for child in placeholder_feature.children:
            #it's possible that the child is also a placeholder, whose attributes haven't been calculated yet!
            if child.start is None:
                self._fixPlaceholder(child)
            child_start_positions.append(child.start)
            child_end_positions.append(child.end)
            
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
        
   
           
    def _mergeCDSSequences(self):
        """GFF files write CDS sequences as separate entries. For DDBJ, we will need to merge them into one single entry,
        so that we can then get the location as join(x..y,X..Y)"""
        
        keys_to_remove = set()
        entries_to_add = []
        
        for key, feature in self.features.items():
            if feature.gfftype == "mRNA":
                cds_list = feature.getAllDownstreamCDS()
                
                if len(cds_list)>1:
                    compound_feature = CompoundFeature(cds_list)
                    entries_to_add.append((compound_feature.getAttribute("ID"), compound_feature))
                    for cds in cds_list:
                        feature.removeDownstreamChild(cds)
                        keys_to_remove.add(cds.attributes["ID"])
                    feature.children.append(compound_feature)
    
        for ktr in keys_to_remove:
            self.features.pop(ktr)
            
        for eta in entries_to_add:
            self.features[eta[0]] = eta[1]
    
    
    
    
    def _detectIncompleteCDS(self):
        """In most GFF files, the start and stop codons are equal to the first and last codon of a CDS.
        However in Breaker2 files, the special features start_codon and stop_codon are provided.
        If i.e. the start_codon is missing, then this means that the CDS sequence is incomplete and this
        needs to be annotated differently in DDBJ annotations.
        """
        if not Parameters.gff_contains_startcodons:
            return
        
        to_replace=[]
        
        for key, feature in self.features.items():
            has_startcodon = False
            has_stopcodon = False
            
                
            if feature.gfftype == 'CDS':
                cds = feature
                parent = feature.parent
                if parent is None:
                    continue
                
                startcodons = parent.getAllDownstreamOfType("start_codon")
                if len(startcodons) > 0:
                    has_startcodon = True
                
                stopcodons = parent.getAllDownstreamOfType("stop_codon")
                if len(stopcodons) > 0:
                    has_stopcodon = True
                
                if len(stopcodons)>1 or len(startcodons)>1:
                    print("ERROR: Multiple start- or stopcodons found for CDS")
                
                if has_startcodon and has_stopcodon: #the CDS is fine, nothing to do here
                    continue
                
                elif has_startcodon: #the end of the CDS is missing

                    if cds.strand == "-":
                        if isinstance(feature, CompoundFeature):
                            #for compound features, we can simply replace the member inquestion
                            cds.members[0] = TruncatedLeftFeature.cloneFeature(cds.members[0])
                            cds._calculatePhase()
                        else:
                            #for regular features, we need to replace them by overwriting the entry in the features
                            #dict later
                            newfeature = TruncatedLeftFeature.cloneFeature(cds)
                            to_replace.append((key, newfeature)) 
                            
                    else: #+strand
                        if isinstance(feature, CompoundFeature):
                            cds.members[-1] = TruncatedRightFeature.cloneFeature(cds.members[-1])
                            cds._calculatePhase()
                        else:
                            newfeature = TruncatedRightFeature.cloneFeature(cds)
                            to_replace.append((key, newfeature))
                    
                        
                elif has_stopcodon:
                    
                    if cds.strand == "-":
                        if isinstance(feature, CompoundFeature):
                            cds.members[-1] = TruncatedRightFeature.cloneFeature(cds.members[-1])
                            cds._calculatePhase()
                        else:
                            newfeature = TruncatedRightFeature.cloneFeature(cds)
                            to_replace.append((key, newfeature))
                    else:
                        if isinstance(feature, CompoundFeature):
                            cds.members[0] = TruncatedLeftFeature.cloneFeature(cds.members[0])
                            cds._calculatePhase()
                        else:
                            newfeature = TruncatedLeftFeature.cloneFeature(cds)
                            to_replace.append((key, newfeature))
                    
                    
                else: #the feature lacks both start and stopcodon
                    if isinstance(feature, CompoundFeature):
                        cds.members[0] = TruncatedLeftFeature.cloneFeature(cds.members[0])
                        cds.members[-1] = TruncatedRightFeature.cloneFeature(cds.members[-1])
                    else:
                        newfeature = TruncatedBothSidesFeature.cloneFeature(cds)
                        to_replace.append((key, newfeature))
                
        for key, newfeature in to_replace:
            oldfeature = self.features[key]
            if oldfeature.parent is not None:
                oldfeature.parent.children.remove(oldfeature)
                oldfeature.parent.children.append(newfeature)
                
            self.features[key] = newfeature
            
                
                
                
                