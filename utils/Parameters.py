from utils.UserInputQuery import UserInputQuery

class Parameters:
    
    @staticmethod
    def init():
        Parameters.params = dict()
        Parameters.qualifiers = set()
        Parameters.string = ""
        Parameters.source_attributes = {"organism":None, "mol_type":None}
        Parameters.keywords = []
    
    @staticmethod
    def addCommonParam(feature_col, qualifier, value):
        feature_dict = Parameters.params.get(feature_col)
        if feature_dict is None:
            feature_dict = dict()
            Parameters.params[feature_col] = feature_dict
        
        values = feature_dict.get(qualifier)
        if values is None:
            values = []
            feature_dict[qualifier] = values
        if value not in values:
            values.append(value)
    
    @staticmethod
    def addCommonParamIfMissing(feature_col, qualifier, value):
        if not Parameters.hasCommonParam(feature_col, qualifier, value):
            Parameters.addCommonParam(feature_col, qualifier, value)
    
    @staticmethod
    def getCommonParams(feature_col, qualifier):
        try:
            return Parameters.params[feature_col][qualifier]
        except:
            return None
    
    @staticmethod
    def hasCommonParam(feature_col, qualifier, value):    
        try:
            values = Parameters.params[feature_col][qualifier]
            return value in values
        except:
            return False
        
        
    @staticmethod
    def parseHeaderFile(header_file):
        try:
            Parameters.params
        except:
            Parameters.init()
        
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
                
                Parameters.addCommonParam(current_feature, qualifier, value)
                
                if qualifier == 'organism':
                    Parameters.source_attributes['organism'] = value
                elif qualifier == 'mol_type':
                    Parameters.source_attributes['mol_type'] = value
        
        
        if Parameters.hasCommonParam("DATATYPE", "type", 'WGS') or Parameters.hasCommonParam("DATATYPE", "type", 'TPA-WGS'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "WGS")
            Parameters.source_attributes["ff_definition"] = "@@[organism]@@ @@[strain]@@ DNA, @@[submitter_seqid]@@"
            
            
        if Parameters.hasCommonParam("DATATYPE", "type", 'TPA') or Parameters.hasCommonParam("DATATYPE", "type", 'TPA-WGS'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "TPA, Third Party Data")
        
        if Parameters.hasCommonParam("DATATYPE", "type", 'TLS'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "TLS, Targeted Locus Study")
        
        if Parameters.hasCommonParam("DIVISION", "division", 'ENV'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "ENV")
        
        if Parameters.hasCommonParam("DIVISION", "division", 'EST'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "EST")
            Parameters.source_attributes["ff_definition"] = "@@[organism]@@ "+UserInputQuery().askForEST2()+", clone: @@[clone]@@"
        
        if Parameters.hasCommonParam("DIVISION", "division", 'GSS'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "GSS")
            Parameters.source_attributes["ff_definition"] = "@@[organism]@@ DNA, clone: @@[clone]@@"
        
        if Parameters.hasCommonParam("DIVISION", "division", 'HTC'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "HTC")
        
        if Parameters.hasCommonParam("DIVISION", "division", 'HTG'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "HTG")
            Parameters.source_attributes["ff_definition"] = "@@[organism]@@ DNA, chromosome @@[map]@@, [BAC/YAC] clone: @@[clone]@@, *** SEQUENCING IN PROGRESS ***"
        
        if Parameters.hasCommonParam("DIVISION", "division", 'STS'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "STS")
            Parameters.source_attributes["ff_definition"] = "@@[organism]@@ DNA, @@[map]@@"
        
        if Parameters.hasCommonParam("DIVISION", "division", 'TSA'):
            Parameters.addCommonParamIfMissing("KEYWORD", "keyword", "TSA, Transcriptome Shotgun Assembly")
        
        
    @staticmethod
    def printCommonParameters():
        for feature_name in Parameters.params:
            print()
            print(feature_name)
            for qualifier in Parameters.params[feature_name]:
                values = Parameters.getCommonParams(feature_name, qualifier)
                for value in values:
                    print("\t"+qualifier, value)
    
    
    
    @staticmethod
    def isInHeader(string):
        try:
            return string in Parameters.string
        except:
            return False
    
    @staticmethod
    def hasQualifier(q):
        return q in Parameters.qualifiers
    
    
    @staticmethod
    def askUserForRequiredParameters():
        userinputquery = UserInputQuery()
        
        #Check if organism and mol_type were present in the COMMON section, since they are required in every 'source' entry
        if not Parameters.hasQualifier('organism'):
            Parameters.source_attributes["organism"] = userinputquery.askUserForOrganism()
        if not Parameters.hasQualifier('mol_type'):
            Parameters.source_attributes["mol_type"] = userinputquery.askUserForMolType()
        
        if Parameters.hasCommonParam("KEYWORD", "keyword", "WGS"):
            if not Parameters.hasCommonParam("KEYWORD", 'keyword', "STANDARD_DRAFT"):
                if not Parameters.hasCommonParam("KEYWORD", 'keyword', "HIGH_QUALITY_DRAFT"):
                    if not Parameters.hasCommonParam("KEYWORD", 'keyword', "IMPROVED_HIGH_QUALITY_DRAFT"):
                        if not Parameters.hasCommonParam("KEYWORD", 'keyword', "NON_CONTIGUOUS_FINISHED"):
                            Parameters.addCommonParam("KEYWORD", "keyword", userinputquery.askForWGS())
            
            if Parameters.source_attributes.get("strain") is None:
                Parameters.source_attributes["strain"] = userinputquery.askUserForStrain()
            
        if Parameters.hasCommonParam("KEYWORD", "keyword", "TPA, Third Party Data"):
            if not Parameters.hasCommonParam("KEYWORD", "keyword", "TPA:inferential"):
                if not Parameters.hasCommonParam("KEYWORD", "keyword", "TPA:experimental"):
                    Parameters.addCommonParam("KEYWORD", "keyword", userinputquery.askForTPA())
            
        if Parameters.hasCommonParam("KEYWORD", 'keyword', 'EST'):
            if not Parameters.hasCommonParam("KEYWORD", "keyword", "5’-end sequence (5’-EST)"):
                if not Parameters.hasCommonParam("KEYWORD", "keyword", "3’-end sequence (3’-EST)"):
                    Parameters.addCommonParam("KEYWORD", "keyword", userinputquery.askForEST())
            
            
            
            
            