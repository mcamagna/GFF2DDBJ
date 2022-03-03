
class Parameters:
    
    
    @staticmethod
    def parseHeaderFile(header_file):
        try:
            Parameters.params
        except:
            Parameters.params = dict()
            Parameters.qualifiers = set()
            Parameters.string = ""
            Parameters.source_attributes = {"organism":None, "mol_type":None}
        
        with open(header_file, 'rt') as filehandle:
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
                    Parameters.string+= current_feature+' '
                qualifier = spl[3]
                value = spl[4]
                
                Parameters.qualifiers.add(qualifier)
                Parameters.string += qualifier +' '+value +' '
                
                current_dict = Parameters.params.get(current_feature)
                if current_dict is None:
                    current_dict = dict()
                    Parameters.params[current_feature] = current_dict
                current_dict[qualifier] = value
                Parameters.params[current_feature] = current_dict
                
                if qualifier == 'organism':
                    Parameters.source_attributes['organism'] = value
                elif qualifier == 'mol_type':
                    Parameters.source_attributes['mol_type'] = value
             
        
    @staticmethod
    def printCommonParameters():
        for feature_name in Parameters.params:
            print()
            print(feature_name)
            for qualifier in Parameters.params[feature_name]:
                print("\t"+qualifier, Parameters.params[feature_name][qualifier])
    
    
    
    @staticmethod
    def isInHeader(string):
        try:
            return string in Parameters.string
        except:
            return False
    
    @staticmethod
    def hasQualifier(q):
        return q in Parameters.qualifiers
    