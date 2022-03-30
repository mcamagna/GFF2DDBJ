'''
@author: Maurizio Camagna
'''
from utils.Parameters import Parameters

def writeGFF(feature_dict: dict):
    """Writes the features from a feature dicht into a GFF file"""
    #sort the dictionary
    features = list(feature_dict.values())
    features.sort(key=lambda x: x.start)
    with open(Parameters.intermediate_gff, 'wt') as out:
        for f in features:
            s = f"{f.seqid}\t{f.source}\t{f.gfftype}\t{f.start}\t{f.end}\t{f.score}\t{f.strand}\t{f.phase}\t"
            for key, value in f.attributes.items():
                s+=f"{key}={value};"
            if s.endswith(";"):
                s = s[:-1]
            s+= "\n"
            
            out.write(s)
      
    
