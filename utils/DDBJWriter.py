'''
@author: Maurizio Camagna
'''
from utils.Parameters import Parameters


class DDBJWriter:
    
    def __init__(self, outpath):
        self.outpath = outpath

        
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
                
                
                current_dict = Parameters.params.get(current_feature)
                if current_dict is None:
                    current_dict = dict()
                    Parameters.params[current_feature] = current_dict
                current_dict[qualifier] = value
                Parameters.params[current_feature] = current_dict
                     
                    
    def writeHeader(self):
        s = "COMMON"

        for feature_name in Parameters.params:
            s+= "\t"
            s+= feature_name
            s+= '\t\t'
            
            for qualifier in Parameters.params[feature_name]:
                values = Parameters.getCommonParams(feature_name, qualifier)
                for value in values:
                    s += qualifier
                    s += '\t'
                    s += value
                    s += '\n\t\t\t'
            #remove last two tabs
            s = s[0:-3]
        with open(self.outpath, 'wt') as out:
            out.write(s)
        
    
    def buildLocationString(self, f):
        """Uses the start/end/strand values to build a location string"""
        s = str(f.start) + '..'+str(f.end)
        if f.strand == '-':
            s= "complement("+s+")"
        return s
    
    
    def _writeFeature(self, f, isSourceFeature=False):
        s = ""
        if isSourceFeature:
            s+= f.seqid + '\t'
            if f.attributes.get("organism") is None:
                f.attributes['organism'] = Parameters.source_attributes['organism']
            if f.attributes.get("mol_type") is None:
                f.attributes['mol_type'] = Parameters.source_attributes['mol_type']
                
        else:
            s += '\t'
        
        s += f.gfftype
        s += '\t'
        s += self.buildLocationString(f)
        s += '\t'
        i = 0
        
        
        for qualifier in f.attributes:
            if i>0:
                s += '\t\t\t'
            value = f.attributes[qualifier]
            s += qualifier + '\t' + value +'\n'
            i+=1
            
        if s[-1] != '\n':
            s+="\n"
            
        with open(self.outpath, 'at') as out:
            out.write(s)
        
        
    def writeFeatures(self, features_dict, sorted_source_feature_keys):
        """Writes a all feature to file. The sorted source features must be provided.
        Note: DDBJ appears to insist that the order of contigs/chromosomes must be the same as 
        in the corresponding fasta file. """
        #source_feature_keys = []
        #for fk in features_dict.keys():
        #    if features_dict[fk].gfftype == "source":
        #        source_feature_keys.append(fk)
        source_feature_keys = sorted_source_feature_keys
        
        for sk in source_feature_keys:
            source_feature = features_dict.get(sk)
            if source_feature is None:
                #print(f"Warning: {sk} does not have a source feature.")
                continue
            if Parameters.sort_features:
                source_feature.sortChildrenByPosition()
                
            self._writeFeature(source_feature, isSourceFeature=True)
            for child in source_feature.children:
                self._writeFeature(child)
        