'''
Braker2 generated GFF files lack the gene feature. We'll try to parse the file as generic as possible to allow conversion of other GFF files as well
@author: Maurizio Camagna
'''
from utils.DDBJWriter import DDBJWriter
from utils.Feature import Feature
from utils.Parameters import Parameters
from utils.CompoundFeature import CompoundFeature

    
    
class GFFParser:
    
    def __init__(self, gff_path):
        self.gff_path = gff_path
        self.features = {} #contains all features
        self.parentFeatures = [] #contains only features that have no parent themselves
    
        self._parseGFF()
        self._connectFeatures()
        self._mergeCDSSequences()
    
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
            
            if not Parameters.isBreaker2_file and source.upper() == "AUGUSTUS":
                Parameters.isBreaker2_file = True
            
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
        
        placeholders = [] #a list of placeholder features, which we may need if some parent features are missing
        
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
        """Placeholder features are features, that were created because the entries in the GFF
        file were referring to a parent, that was not written into the GFF file.
        Using their child nodes, we can infer the required information.
        """
        children_contain_gene_feature = False
        child_gene = None
        
        
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
                child_gene = child
            
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
            #TODO: Actually the placeholder is probably a transcript, not a gene.
            #We can treat it as a gene, since these transcripts will later be merged anyways,
            #but for compatibility with non-braker2 gff files, this needs to be fixed
            placeholder_feature.gfftype = "gene"
        else:
            #The placeholder contains a gene feature as child.
            #We will make that gene the parent of the other children of the placeholder
            #before removing the placeholder, since it serves no purpose
            new_parent = child_gene
            for pc in placeholder_feature.children:
                if pc != new_parent:
                    new_parent.children.add(pc)
                    pc.parent = new_parent
            self.parentFeatures.remove(placeholder_feature)
            self.parentFeatures.append(new_parent)
            self.features.pop(placeholder_feature.attributes.get("ID"))
        
           
    def _mergeCDSSequences(self):
        """GFF files write CDS sequences as separate entries. For DDBJ, we will need to merge them into one single entry,
        so that we can then get the location as join(x..y,X..Y)"""
        
        keys_to_remove = set()
        entries_to_add = []
        
        for key, feature in self.features.items():
            if feature.gfftype == "gene":
                cds_list = feature.getAllDownstreamCDS()
                
                if len(cds_list)>0:
                    compound_feature = CompoundFeature(cds_list)
                    entries_to_add.append((compound_feature.getAttribute("ID"), compound_feature))
                    for cds in cds_list:
                        feature.children.remove(cds)
                        keys_to_remove.add(cds.attributes["ID"])
                    feature.children.append(compound_feature)
    
        for ktr in keys_to_remove:
            self.features.pop(ktr)
            
        for eta in entries_to_add:
            self.features[eta[0]] = eta[1]
    
    
    