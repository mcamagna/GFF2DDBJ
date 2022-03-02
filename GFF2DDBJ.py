'''
@author: Maurizio Camagna
'''
from utils.GFFParser import GFFParser



parser = GFFParser("data/Aiko/Fusarium_genome/braker_217701/braker.gff3")
print("Features:", len(parser.features))




OUTFILE = "test.txt"





#parser = GFFParser("/mnt/Data/GenomeData/Medicago/Annotation/Mt4.0v2_genes_20140818_1100.gff3")
#print("Features:", len(parser.features))

#parser = GFFParser("/mnt/Data/GenomeData/Nicotiana_benthamiana/annotation/Niben101_annotation.allfeatures.gff.gz")
#print("Features:", len(parser.features))