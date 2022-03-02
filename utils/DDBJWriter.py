'''
@author: Maurizio Camagna
'''

class DDBJWriter:
    
    def __init__(self, outpath, header_file):
        self.outpath = outpath
        self.header_file = header_file
        self.common = dict()
        
        if self.header_file is not None:
            self._parseHeaderFile()
        
        self.source_attributes = {"organism":None, "mol_type":None}
        #search the parsed header whether they contain the mandatory values for organism and mol_type
        for _dict in self.common.values():
            if "organism" in _dict.keys():
                self.source_attributes["organism"] = _dict.get("organism")
            if "mol_type" in _dict.keys():
                self.source_attributes["mol_type"] = _dict.get("mol_type")

        
    def _parseHeaderFile(self):
        with open(self.header_file, 'rt') as filehandle:
            current_feature = ""
            
            for line in filehandle:
                line = line.replace("\n", '')
                #skip comments
                if line.startswith("#"):
                    continue
                #skip empty lines
                if len(line.replace(" ", "").replace("\t", "")) ==0:
                    continue
            
                spl = line.split("\t")
                if spl[1] != "":
                    current_feature = spl[1]
                qualifier = spl[3]
                value = spl[4]
                
                
                current_dict = self.common.get(current_feature)
                if current_dict is None:
                    current_dict = dict()
                    self.common[current_feature] = current_dict
                current_dict[qualifier] = value
                self.common[current_feature] = current_dict
            
    def organismWasProvided(self):
        return self.source_attributes.get("organism") is not None
                    
    def molTypeWasProvided(self):
        return self.source_attributes.get("mol_type") is not None
                    
                    
    def writeHeader(self):
        s = "COMMON"

        for feature_name in self.common:
            s+= "\t"
            s+= feature_name
            s+= '\t\t'
            
            for qualifier in self.common[feature_name]:
                value = self.common[feature_name][qualifier]
                s += qualifier
                s += '\t'
                s += value
                s += '\n\t\t\t'
            #remove last two tabs
            s = s[0:-3]
        with open(self.outpath, 'wt') as out:
            out.write(s)
        
        
    
    def printCommonParameters(self):
        for feature_name in self.common:
            print()
            print(feature_name)
            for qualifier in self.common[feature_name]:
                print("\t"+qualifier, self.common[feature_name][qualifier])
    
    
    
    def buildLocationString(self, f):
        """Uses the start/end/strand values to build a location string"""
        s = str(f.start) + '..'+str(f.end)
        if f.strand == '-':
            s= "complement("+s+")"
        return s
    
    
    def writeFeature(self, f, isSourceFeature=False):
        s = ""
        if isSourceFeature:
            s+= f.seqid + '\t'
            if f.attributes.get("organism") is None:
                f.attributes['organism'] = self.source_attributes['organism']
            if f.attributes.get("mol_type") is None:
                f.attributes['mol_type'] = self.source_attributes['mol_type']
        else:
            s += '\t'
        
        s += f.gfftype
        s += '\t'
        s += self.buildLocationString(f)
        s += '\t'
        i = 0
        
        #TODO: All GFF attributes need to be checked whether they are permissible in a DDBJ annotation
        attr = f.attributes
        if attr.get("ID") is not None:
            attr.pop("ID")
        if attr.get("Parent") is not None:
            attr.pop("Parent")
        ###############################
        
        for qualifier in attr:
            if i>0:
                s += '\t\t\t'
            value = f.attributes[qualifier]
            s += qualifier + '\t' + value +'\n'
            i+=1
        
        if s[-1] != '\n':
            s+="\n"
            
        with open(self.outpath, 'at') as out:
            out.write(s)
        
        
         
    def writeWholeFeature(self, feature):
        """Writes a feature, including subfeatures such as CDS, UTR's etc"""
        if feature.gfftype.lower() != "source":
            source_feature = feature.clone()
            source_feature.gfftype = 'source'
            self.writeFeature(source_feature, isSourceFeature=True)
        self.writeFeature(feature)
        
        for child in feature.children:
            child_type = child.gfftype.lower()
            if child_type == 'cds' or child_type == "exon" or child_type == 'mrna':
                self.writeFeature(child)
        pass