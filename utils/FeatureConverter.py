'''
This will help to check whether GFF features are permissable as DDBJ entries and attempt to modify them appropriately.
@author: Maurizio Camagna
'''
from utils.GFFParser import Feature
from utils.DDBJWriter import DDBJWriter

class FeatureConverter:
    
    def __init__(self):
        self.parseFeatureList()
        #in order to easily find features that have a mismatch in lower/upper case
        #we'll also make a dict that maps alternative names to each DDBJ feature
        self.generateFeatureMappings()
        self.generateQualifierMappings()
              
    def parseFeatureList(self):
        """Parses the 'DDBJ_Features.tsv' file, which contains the information
        on allowed DDBJ features and qualifiers.
        """
        self.ddbj_features = dict()
        
        with open("DDBJ_Features.tsv", 'rt') as f:
            for i, line in enumerate(f):
                if i==0:
                    continue
                line = line.replace("\n", "")
                s = line.split("\t")
                feature_dict = self.ddbj_features.get(s[0])
                if feature_dict is None:
                    feature_dict = dict()
                    feature_dict["Mandatory"] = set()
                    feature_dict["Optional"] = set()
                    self.ddbj_features[s[0]] = feature_dict
                
                if s[1] == 'M':
                    mandatory_set = feature_dict.get("Mandatory")
                    mandatory_set.add(s[2])
                    
                elif s[1] == 'O':
                    optional_set = feature_dict.get("Optional")
                    optional_set.add(s[2])
    
    
    def generateFeatureMappings(self):
        """Generates a dict that maps various spelling variations to the correct DDBJ Feature type"""
        self.ddbj_feature_mappings = dict()
        
        self.ddbj_feature_mappings["cDNA"] = 'mRNA'
        self.ddbj_feature_mappings["cdna"] = 'mRNA'
        self.ddbj_feature_mappings["CDNA"] = 'mRNA'
        
        for k in self.ddbj_features.keys():
            alt_names = set()
            alt_names.add(k)
            alt_names.add(k.replace("'", "-"))
            alt_names.add(k.replace("-", ""))
            alt_names.add(k.replace("-", "_"))
            alt_names.add(k.replace("_", "-"))
            alt_names.add(k.replace("_", ""))
            alt_names.add(k.replace("3'UTR", "three_prime_UTR"))
            alt_names.add(k.replace("3", "three"))
            alt_names.add(k.replace("5'UTR", "five_prime_UTR"))
            alt_names.add(k.replace("5", "five"))
            
            for alt_name in alt_names:
                self.ddbj_feature_mappings[alt_name] = k
            for alt_name in alt_names:
                self.ddbj_feature_mappings[alt_name.lower()] = k  
    
    
    
    def generateQualifierMappings(self):
        """Generates a dict that maps various spelling variations to the correct DDBJ qualifier"""
        self.ddbj_qualifier_mappings = dict()
        qualifiers = set()
        for _dict in self.ddbj_features.values():
            qualifiers.update(_dict["Mandatory"])
            qualifiers.update(_dict["Optional"])
            
        for qualifier in qualifiers:
            self.ddbj_qualifier_mappings[qualifier] = qualifier
            self.ddbj_qualifier_mappings[qualifier.replace("_", "-")] = qualifier
            self.ddbj_qualifier_mappings[qualifier.replace("_", "")] = qualifier
            self.ddbj_qualifier_mappings[qualifier.lower()] = qualifier
            self.ddbj_qualifier_mappings[qualifier.lower().replace("_", "-")] = qualifier
            self.ddbj_qualifier_mappings[qualifier.lower().replace("_", "")] = qualifier
        
        self.ddbj_qualifier_mappings["ID"] = "locus_tag"
        self.ddbj_qualifier_mappings["id"] = "locus_tag"
        self.ddbj_qualifier_mappings["Id"] = "locus_tag"
        
    
    def _addSourceFeatures(self, gff_feature_dict):
        """DDBJ annotation files require a 'source' feature. Typically, this source represents the chromosome or contig
        and we can therefore use the GFF seqid to merge GFF entries that belong to the same 'source'
        """
        grouped_locations = dict()
        
        for key, feature in gff_feature_dict.items():
            group = grouped_locations.get(feature.seqid)
            if group is None:
                group = {"name":feature.seqid, "features":[], "start_list":[], "end_list":[]}
                grouped_locations[feature.seqid] = group
                
            group["features"].append(feature)
            group["start_list"].append(feature.start)
            group["end_list"].append(feature.end)
        
        for group_key, group in grouped_locations.items():
            #start = min(group["start_list"])
            #UME validation showed an error. It seems like the source start must always equal 1?!
            start = 1
            
            end = max(group["end_list"])
            attr = DDBJWriter.source_attributes.copy()
            group_feature = Feature(seqid=group_key,gfftype="source", start=start, end=end, attribute_dict=attr)
            
            #Drop all hierarchies in all subnodes and make the source feature the parent of all children
            group_feature.children = group["features"]
            for child in group_feature.children:
                child.parent = group_feature
                child.children.clear()
            gff_feature_dict[group_key] = group_feature
     
     
    def _removeGeneFeatures(self, gff_feature_dict):
        """Genes are not allowed in DDBJ annotation files. Instead of simply removing them,
        we should try and pass the information of the Gene on to it's child features if possible."""
        to_remove_keys = set()
        for key, feature in gff_feature_dict.items():
            if feature.gfftype.lower() == "gene":
                to_remove_keys.add(key)
                #let's transfer the gene annotation to all downstream CDS and mRNA's
                child_features_to_augment = set()
                child_features_to_augment.update(feature.getAllDownstreamCDS())
                child_features_to_augment.update(feature.getAllDownstreamMRNA())
                
                for child in child_features_to_augment:
                    #copy all gene attributes from to CDS/mRNA if no conflicting attribute is present
                    for gene_attr in feature.attributes.keys():
                        if gene_attr not in child.attributes.keys():
                            child.attributes[gene_attr] = feature.attributes[gene_attr]
                    #add the 'gene' qualifier to the child.
                    if "gene" not in child.attributes.keys():
                        if "standard_name" in feature.attributes.keys():
                            child.attributes["gene"] = feature.attributes["standard_name"]
                        elif "locus_tag" in feature.attributes.keys():
                            child.attributes["gene"] = feature.attributes["locus_tag"]
                            
                #let's dissolve the child/parent relationships for the gene node
                for child in feature.children:
                    child.parent = feature.parent
                feature.children = None #unnecessary but throws an error if the feature is accidentally used elsewhere
        #Let's remove the gene feature
        [gff_feature_dict.pop(r) for r in to_remove_keys]
            
    
    def convertFeatures(self, gff_feature_dict):
        len_before = len(gff_feature_dict)
        
        self._mapGFF_Features(gff_feature_dict)
        self._mapQualifiers(gff_feature_dict)
        self._removeGeneFeatures(gff_feature_dict)
        self._addSourceFeatures(gff_feature_dict)
        self._checkValidityOfQualifiers(gff_feature_dict)
        self._removeEntriesWithouthQualifiers(gff_feature_dict)
        
        
        len_after = len(gff_feature_dict)
        removed_feature_count = len_before-len_after
        if removed_feature_count>0:
            print(f"Number of invalid GFF entries that will not be converted: {removed_feature_count}")
    
    def _removeEntriesWithouthQualifiers(self, gff_feature_dict):
        to_remove = set()
        
        for fkey in gff_feature_dict.keys():
            attr = gff_feature_dict[fkey].attributes
            if len(attr)==0:
                to_remove.add(fkey)

        [gff_feature_dict.pop(r) for r in to_remove]
        
    def _checkValidityOfQualifiers(self, gff_feature_dict):
        for fkey in gff_feature_dict.keys():
            gff_feature = gff_feature_dict[fkey]
            filtered_attributes = dict()
            for qualifier in gff_feature.attributes.keys():
                if qualifier in self.ddbj_features[gff_feature.gfftype]["Mandatory"] or qualifier in self.ddbj_features[gff_feature.gfftype]["Optional"]:
                    filtered_attributes[qualifier] = gff_feature.attributes[qualifier]
            gff_feature.attributes = filtered_attributes
            
            all_mandatory_present = True
            for mandatory_qualifier in self.ddbj_features[gff_feature.gfftype]["Mandatory"]:
                if mandatory_qualifier not in filtered_attributes.keys():
                    all_mandatory_present = False
            if not all_mandatory_present:
                print(f"ERROR: Mandatory qualifier missing in GFF type {gff_feature}.\nDDBJ requires the following qualifiers for this feature:", self.ddbj_features[gff_feature.gfftype]["Mandatory"])
                import sys
                sys.exit(1)
                
                
    def _mapQualifiers(self, gff_feature_dict):
        invalid_qualifiers = set()
        
        for fkey in gff_feature_dict.keys():
            gff_feature = gff_feature_dict[fkey]
            converted_attributes = dict()
            
            for qualifier in gff_feature.attributes.keys():
                converted_qualifier = self._convertQualifier(qualifier)
                if converted_qualifier is None:
                    invalid_qualifiers.add(qualifier)
                else:
                    converted_attributes[converted_qualifier] = gff_feature.attributes[qualifier]
            
            gff_feature.attributes = converted_attributes
            
            
    def _convertQualifier(self, gff_qualifier):
        """"""
        hit = self.ddbj_qualifier_mappings.get(gff_qualifier)
        if hit is None:
            hit = self.ddbj_qualifier_mappings.get(gff_qualifier.lower())
            if hit is None:
                return hit
        return hit
                
            
    
    def _mapGFF_Features(self, gff_feature_dict):
        """Checks GFF Feature types and removes/renames invalid feature types"""
        invalid_gff_feature_types = set()
        gff_features_to_remove = set()
        
        for fkey in gff_feature_dict.keys():
            gff_feature = gff_feature_dict[fkey]
            
            converted_type = self._convertGFF_FeatureType(gff_feature.gfftype)
            if converted_type is None:
                gff_features_to_remove.add(fkey)
                invalid_gff_feature_types.add(gff_feature.gfftype) 
                continue
            else:
                gff_feature.gfftype = converted_type
        
        #remove GFF entries with invalid types
        [gff_feature_dict.pop(r) for r in gff_features_to_remove]
            
            #attr_dict = gff_feature.attributes
            #for akey in gff_feature.attributes.keys():
        if len(invalid_gff_feature_types)>0:
            print("WARNING: The following invalid feature types will be omitted: ", invalid_gff_feature_types)
    
    
    
    def _convertGFF_FeatureType(self, gfftype):
        hit = self.ddbj_feature_mappings.get(gfftype)
        if hit is None:
            hit = self.ddbj_feature_mappings.get(gfftype.lower())
            return hit
        return hit   
    
    
    
    
    