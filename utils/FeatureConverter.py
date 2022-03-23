'''
This will help to check whether GFF features are permissable as DDBJ entries and attempt to modify them appropriately.
@author: Maurizio Camagna
'''

from utils.FastaParser import FastaParser
from utils.Parameters import Parameters
from utils.features import Feature, CompoundFeature, TruncatedStartFeature, TruncatedEndFeature,TruncatedFeature, TruncatedBothSidesFeature
import re

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
        
        #self.ddbj_qualifier_mappings["ID"] = "locus_tag"
        #self.ddbj_qualifier_mappings["id"] = "locus_tag"
        #self.ddbj_qualifier_mappings["Id"] = "locus_tag"
        
    
    def _addSourceFeatures(self, gff_feature_dict):
        """DDBJ annotation files require a 'source' feature. Typically, this source represents the chromosome or contig
        and we can therefore use the GFF seqid to obtain the data from the fasta file.
        """
        needs_submitter_seqid = False
        try:
            datatype = Parameters.getCommonParams("DATATYPE", "type")
            if "WGS" in datatype or "TSA" in datatype or "TLS" in datatype or "CON" in datatype:
                needs_submitter_seqid = True
        except:
            pass 
        
        grouped_locations = dict()
        
        for key, feature in gff_feature_dict.items():
            group = grouped_locations.get(feature.seqid)
            if group is None:
                group = {"name":feature.seqid, "features":[]}
                grouped_locations[feature.seqid] = group
            group["features"].append(feature)
        
        for group_key, group in grouped_locations.items():
            start = 1
            end = FastaParser.fasta_dict[group_key] -2 #-2 because this is what UME validator asks for
            attr = Parameters.source_attributes.copy()
            if needs_submitter_seqid:
                attr["submitter_seqid"] = group_key
                
            group_feature = Feature(seqid=group_key,gfftype="source", start=start, end=end, attribute_dict=attr)
            
            #Drop all hierarchies in all subnodes and make the source feature the parent of all children
            group_feature.children = group["features"]
            for child in group_feature.children:
                child.parent = group_feature
                child.children.clear()
            gff_feature_dict[group_key] = group_feature
        
        #finally we need to also add all contigs/chromosomes that
        #are present in the fasta file, but not in the GFF file
        remaining = set(FastaParser.fasta_dict.keys()).difference(set(grouped_locations.keys()))
        for key in remaining:
            start = 1
            end = FastaParser.fasta_dict[key] -2 #-2 because this is what UME validator asks for (lenght starts at 1 and includes last)
            attr = Parameters.source_attributes.copy()
            if needs_submitter_seqid:
                attr["submitter_seqid"] = key
            feature = Feature(seqid=key,gfftype="source", start=start, end=end, attribute_dict=attr)
            gff_feature_dict[key] = feature
    
    
    def _removePlaceHolderTranscriptsFeatures(self, gff_feature_dict):
        """During the parsing of the GFF file, we may have added mRNA features that were never present in the GFF file.
        We'll now remove these mRNA's, since writing them to the annotation could lead to the wrong conclusion that these
        mRNA's experimentally obtained """
        if Parameters.gff_contains_transcripts:
            #the original gff file provided transcripts. We shouldn't remove the mRNAs
            return
        
        to_remove_keys = set()
        for key, feature in gff_feature_dict.items():
            if feature.gfftype == "mRNA":
                to_remove_keys.add(key)
                for child in feature.children:
                    child.parent = feature.parent
                    
        [gff_feature_dict.pop(r) for r in to_remove_keys]
     
    def _removeGeneFeatures(self, gff_feature_dict):
        """Genes are not allowed in DDBJ annotation files. Instead of simply removing them,
        we should try and pass the information of the Gene on to it's child features if possible."""
        to_remove_keys = set()
        for key, feature in gff_feature_dict.items():
            if feature.gfftype.lower() == "gene":
                to_remove_keys.add(key)
                #let's transfer the gene annotation to all downstream CDS and mRNA's, introns, exons, etc
                child_features_to_augment = set()
                child_features_to_augment.update(feature.getAllDownstreamChildren())
                for child in child_features_to_augment:
                    #copy all gene attributes from to CDS/mRNA if no conflicting attribute is present
                    for gene_attr in feature.attributes.keys():
                        if gene_attr not in child.attributes.keys():
                            child.attributes[gene_attr] = feature.attributes[gene_attr]
                            
                #let's dissolve the child/parent relationships for the gene node
                for child in feature.children:
                    child.parent = feature.parent
                feature.children = None #unnecessary but throws an error if the feature is accidentally used elsewhere
        #Let's remove the gene feature
        [gff_feature_dict.pop(r) for r in to_remove_keys]
    
    
    def addAssemblyGaps(self, gff_feature_dict, gaps):
        """Adds assembly gaps found in the FASTA file to the feature_dict"""
        for contig_name in gaps.keys():
            gaplist = gaps[contig_name]
            for i, gap in enumerate(gaplist):
                attr = Parameters.assembly_gap_attributes.copy()
                name = contig_name+"_assembly_gap_"+str(i)
                f = Feature(seqid=contig_name, gfftype="assembly_gap", start=gap[0]+1, end=gap[1], strand="+", attribute_dict=attr)
                f.parent = gff_feature_dict.get(contig_name)
                gff_feature_dict.get(contig_name).children.append(f)
                gff_feature_dict[name] = f
                
        self._splitCDSWithGaps(gff_feature_dict)
        #if CDS sequences were split, they will contain the ID attribute again
        #therefore:
        self._checkValidityOfQualifiers(gff_feature_dict)
        self._removeEntriesWithouthQualifiers(gff_feature_dict)
        
        
        
    def _fixLocusTagsAndGeneNames(self, gff_feature_dict):
        """The DDBJ locus tag naming convention is very strict.
        The locus tag must be preceded by a locus tag prefix, separated by an underscore. Locus tags are assigned to most subfeatures 
        of genes and these subfeatures must have the identical locus tag as the corresponding gene BUT the tag cannot be the same as the gene name.
        Also, in case all subfeatures share the same locus_tag and have a gene qualifier, then the locus_tag should be removed in favour of the gene name.
         Note: this function will assign the genes locus tag to all subfeatures, regardless of the type. Invalid assigning of the locus_tag 
         qualifier will need to be filtered out by _checkValidityOfQualifiers()
        """
        prefix = Parameters.locus_attributes["locus_tag_prefix"]
        need_to_add_prefix = prefix!=""
        
        for key, feature in gff_feature_dict.items():
            if feature.gfftype == 'gene':
                locus_tag = feature.attributes.get("locus_tag")
                if locus_tag is None:
                    #need to build a locus tag from the gene name and strip all non-numeric values
                    locus_tag = feature.getAttribute("gene")
                    locus_tag = re.sub('[^0-9]','', locus_tag)
                    #let's pad the number with zeros
                    locus_tag = ("0"*(8-len(locus_tag)))+locus_tag
                    
                
                if need_to_add_prefix:
                    #remove underscores since they are not permissible
                    if "_" in locus_tag:
                        locus_tag = locus_tag.replace("_", "")
                        
                    locus_tag = prefix+"_"+locus_tag
                    feature.attributes["locus_tag"] = locus_tag
                    for child in feature.children:
                        #We will assign both gene and locus tag at this point.
                        #During the _checkValidityOfQualifiers step, the correct choice
                        #the locus tag may be removed depending on the circumstances.
                        child.attributes["gene"] = feature.attributes["gene"]
                        child.attributes["locus_tag"] = locus_tag

    
    def _removeDuplicateFeatures(self, gff_feature_dict):
        feature_keys_to_remove = set()
        hash_set = set()
        
        for key, feature in gff_feature_dict.items():
            feature.removeDuplicateChildren()
            
            h = feature.getHash()
            if h in hash_set:
                feature_keys_to_remove.add(key)
            else:
                hash_set.add(h)
        
        for fkr in feature_keys_to_remove:
            gff_feature_dict.pop(fkr)
    
    
    @staticmethod
    def _calculateORF_Offset(pos_1, pos_2):
        """Calculates the Open Reading Frame offset.
        If pos 1 and pos 2 are in frame, return 0
        otherwise returns 1, or 2.
        The order of the positions is irrelevant
        """
        a = min(pos_1, pos_2)
        b = max(pos_1, pos_2)
        diff = b-a
        diff+=1 #Not sure, but the b position maybe up-to, but not including the last base
        diff+=1 #because we need the number of base pairs. Example 6bp's:  1,2,3,4,5,6 --> 6-1 = 5
        mod = diff%3
        return mod
        
    
    def _addAdditionalCDSQualifiers(self, gff_feature_dict):
        """DDBJ annotations require CDS features to contain transl_table and codon_start.
        codon_start represents the ORF offset and would typically be 1, since a CDS starts with the start codon.
        However, the CDS may also be a truncated CDS with an unknown start or end. We'll try to infer the codon_start
        in such cases.
        transl_table will always be set to 1 if no other value was provided via the command line.
        """
                   
        for key, feature in gff_feature_dict.items():
            if feature.gfftype == "CDS":
                tt = feature.attributes.get("transl_table")
                if tt is None:
                    tt = "1"
                    feature.attributes["transl_table"] = tt
                        
                cs = feature.attributes.get("codon_start")
                if cs is None:
                    
                    cs = str((int(feature.phase)+1))
                    
                    feature.attributes["codon_start"] = cs
    
    def _addExonIntronNumbers(self, gff_feature_dict):
        """Exons and introns need to be numbered by occurence in 5->3 direction."""
        #We will first need to regroup the exons and introns by their corresponding gene tags
        gene_groups = dict()
        for key, feature in gff_feature_dict.items():
            if feature.gfftype == "exon" or feature.gfftype == "intron":
                gene_group = gene_groups.get(feature.getAttribute("gene"))
                if gene_group is None:
                    gene_group = []
                    gene_groups[feature.getAttribute("gene")] = gene_group
                gene_group.append(feature)
                
        for group in gene_groups.values():
            group.sort(key=lambda x: x.start, reverse=True if group[0].strand == "-" else False)
            for i, element in enumerate(group):
                element.attributes["number"] = str(i+1)
               
     
    def _splitCDSWithGaps(self, gff_feature_dict):
        """Finds assembly_gaps that lie within the bounds of a CDS, and splits the CDS into two pieces"""
        
        features_to_add = []
        keys_to_remove = set()
        
        for key, feature in gff_feature_dict.items():
            if feature.gfftype == 'source':
                cds_list = list(feature.getAllDownstreamOfType("CDS"))
                cds_list.sort(key=lambda x: x.start)
                
                gap_list = list(feature.getAllDownstreamOfType("assembly_gap"))
                gap_list.sort(key=lambda x: x.start)
                
                #we will iterate through the sorted CDS and gap list to find overlaps of CDSs and gaps
                gap_index = 0
                cds_index = 0
                while(gap_index<len(gap_list) and cds_index<len(cds_list)):
                    cds = cds_list[cds_index]
                    gap = gap_list[gap_index]
                    if cds.end < gap.start:
                        cds_index+=1
                    elif gap.end < cds.start:
                        gap_index +=1
                    else:
                        #print(f"CDS: {cds.buildLocationString()}\nGAP:{gap.buildLocationString()}")
                        #compound CDS's don't necessarily need a split, since the gap may lie inside an intron
                        split_required = False
                        if isinstance(cds, CompoundFeature):
                            for subcds in cds.members:
                                if subcds.end<gap.start or gap.end<subcds.start:
                                    continue
                                else:
                                    split_required = True
                        else:
                            split_required = True
                            
                        if split_required:
                            #print("Split required!")
                            split_cdss = cds.split(gap.start, gap.end)
                            left = split_cdss[0]
                            right = split_cdss[1]
                            feature.children.remove(cds)
                            feature.children.append(left)
                            feature.children.append(right)
                                
                            #we need the key of the cds feature
                            cds_keys = [k for k, val in gff_feature_dict.items() if val == cds]
                            if len(cds_keys)==0:
                                #the key was not found, because it's a CDS that was previously split and hasn't
                                #been added to the gff_feature_dict yet
                                cds_keys = [k for k, val in features_to_add if val == cds]
                                
                            cds_key = cds_keys[0]
                            keys_to_remove.add(cds_key)
                            features_to_add.append((cds_key+"_l", left))
                            features_to_add.append((cds_key+"_r", right))
                            
                            #also, a single CDS could have multiple gaps. We need to re-generate the cds_list
                            #since the newly split CDS will now be found via getAllDownstreamOfType
                            cds_list = list(feature.getAllDownstreamOfType("CDS"))
                            cds_list.sort(key=lambda x: x.start)
                            gap_index-=1 # we want to rerun the same gap again in case we have multiple CDS's affected by the same gap
                        
                        gap_index+=1
        
        for k in keys_to_remove:
            try:
                gff_feature_dict.pop(k)
            except:
                pass
            
        for key, value in features_to_add:
            gff_feature_dict[key] = value 
          
          
                    
    def convertFeatures(self, gff_feature_dict):
        len_before = len(gff_feature_dict)
        
        self._addAdditionalCDSQualifiers(gff_feature_dict)
        self._mapGFF_Features(gff_feature_dict)
        self._mapQualifiers(gff_feature_dict)
        self._fixLocusTagsAndGeneNames(gff_feature_dict)
        self._removePlaceHolderTranscriptsFeatures(gff_feature_dict)
        self._removeGeneFeatures(gff_feature_dict)
        self._addSourceFeatures(gff_feature_dict)
        self._checkValidityOfQualifiers(gff_feature_dict)
        self._removeEntriesWithouthQualifiers(gff_feature_dict)
        
        #Braker2 was found to annotate the same region multiple times, with slightly different ID's
        #after conversion, it is possible that we end up with identical features. Let's remove them
        self._removeDuplicateFeatures(gff_feature_dict)
        #exon/intron numbers are added after removing of duplicates has succeeded
        #otherwise they would receive different hashes
        self._addExonIntronNumbers(gff_feature_dict)
        
        
        len_after = len(gff_feature_dict)
        removed_feature_count = len_before-len_after
        #if removed_feature_count>0:
        #    print(f"Number of invalid GFF entries that will not be converted: {removed_feature_count}")
        
        return gff_feature_dict
    
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
            
            #If both gene and locus_tag qualifiers are present, keep only the gene qualifier
            if gff_feature.attributes.get("gene") is not None and gff_feature.attributes.get("locus_tag") is not None:
                gff_feature.attributes.pop('locus_tag')
            
            all_mandatory_present = True
            for mandatory_qualifier in self.ddbj_features[gff_feature.gfftype]["Mandatory"]:
                if mandatory_qualifier not in filtered_attributes.keys():
                    all_mandatory_present = False
            if not all_mandatory_present:
                print(f"ERROR: Mandatory qualifier missing in GFF type {gff_feature}.\nDDBJ requires the following qualifiers for this feature:", self.ddbj_features[gff_feature.gfftype]["Mandatory"])
                import sys
                sys.exit(1)
                
                
    def _mapQualifiers(self, gff_feature_dict):
        """Maps/converts GFF qualifiers to DDBJ qualifiers if possible and removes invalid qualifiers otherwise """
        invalid_qualifiers = set()
        
        for fkey in gff_feature_dict.keys():
            gff_feature = gff_feature_dict[fkey]
            converted_attributes = dict()
            
            for qualifier in gff_feature.attributes.keys():
                #Special case: GFF ID's are invalid DDBJ qualifiers, and would be removed,
                #however, we do need to keep this information for genes, since the gene name is required
                #to be passed to child nodes and also for the locus_tag.
                converted_qualifier = None
                if gff_feature.gfftype.lower() == "gene" and qualifier.upper() == "ID" and not gff_feature.hasAttribute("gene"):
                    converted_qualifier = "gene"
                else:
                    converted_qualifier = self._convertQualifier(qualifier)
                
                if converted_qualifier is None:
                    invalid_qualifiers.add(qualifier)
                else:
                    converted_attributes[converted_qualifier] = gff_feature.attributes[qualifier]
            
            gff_feature.attributes = converted_attributes
            
            
    def _convertQualifier(self, gff_qualifier):
        """Looks up GFF qualifier names in the ddbj_qualifier_mappings dict and returns the matching DDBJ qualifier."""
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
    
    
    
    
    