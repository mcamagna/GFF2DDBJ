'''
@author: Maurizio Camagna
'''
import sys
from utils.GFFParser import GFFParser
from utils.DDBJWriter import DDBJWriter
from utils.UserInputQuery import UserInputQuery

#INFILE = "data/Aiko/Fusarium_genome/braker_217701/braker.gff3"
INFILE = "/mnt/Data/GenomeData/Tomato/ITAG4.1_release/ITAG4.1_gene_models.gff"
OUTFILE = "test.txt"
HEADERFILE = "test_header.txt"

ddbjwriter = DDBJWriter(OUTFILE, HEADERFILE)
print("The COMMON header so far contains these values:")
ddbjwriter.printCommonParameters()


userinputquery = UserInputQuery()

#Check if organism and mol_type were present in the COMMON section, since they are required in every 'source' entry
if not ddbjwriter.organismWasProvided():
    ddbjwriter.source_attributes["organism"] = userinputquery.askUserForOrganism()
if not ddbjwriter.molTypeWasProvided():
    ddbjwriter.source_attributes["mol_type"] = userinputquery.askUserForMolType()


print("Parsing GFF file:", INFILE)
gffparser = GFFParser(INFILE)
print("Number of features found in GFF file:", len(gffparser.features))

ddbjwriter.writeHeader()

for parent_feature in gffparser.parentFeatures:    
    ddbjwriter.writeWholeFeature(parent_feature)

print("Conversion finished...")






#parser = GFFParser("/mnt/Data/GenomeData/Medicago/Annotation/Mt4.0v2_genes_20140818_1100.gff3")
#print("Features:", len(parser.features))

#parser = GFFParser("/mnt/Data/GenomeData/Nicotiana_benthamiana/annotation/Niben101_annotation.allfeatures.gff.gz")
#print("Features:", len(parser.features))