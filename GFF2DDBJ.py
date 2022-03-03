'''
@author: Maurizio Camagna
'''
import sys
from utils.GFFParser import GFFParser
from utils.DDBJWriter import DDBJWriter
from utils.UserInputQuery import UserInputQuery
from utils.FeatureConverter import FeatureConverter
from utils.FastaParser import FastaParser

#INFILE = "data/GCF_000143535.2_ASM14353v4_genomic.gff"
#INFILE = "/mnt/Data/GenomeData/Tomato/ITAG4.1_release/ITAG4.1_gene_models.gff"
#INFILE = "/mnt/Data/GenomeData/Eggplant/Eggplant_V4.1_function_IPR_final.gff"
#INFILE = "/mnt/Data/GenomeData/Nicotiana_benthamiana/annotation/Niben101_annotation.gene_models.gff"


#Fusarium langsethiae
INFILE = "data/Aiko/Fusarium_genome/braker_217701/braker.gff3"
FASTAFILE = 'data/Aiko/Fusarium_genome/SAMD00414474_MFG217701.fasta'
OUTFILE = "data/SAMD00414474_MFG217701.ann"
HEADERFILE = "data/aiko_header_217701.txt"

#INFILE = "data/Aiko/Fusarium_genome/braker_217702/braker.gff3"
#FASTAFILE = 'data/Aiko/Fusarium_genome/SAMD00414475_MFG217702.fasta'
#OUTFILE = "data/SAMD00414474_MFG217702.ann"
#HEADERFILE = "data/aiko_header_217702.txt"

#TODO: check if all filepaths are found before running heavy operations

ddbjwriter = DDBJWriter(OUTFILE, HEADERFILE)
print("The COMMON header currently contains these values:")
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
features = gffparser.features

fconverter = FeatureConverter()
fconverter.convertFeatures(features)

ddbjwriter.writeHeader()

fastaParser = FastaParser(FASTAFILE)
fasta_headers = fastaParser.getFastaHeaders()
ddbjwriter.writeFeatures(features, fasta_headers)

print("Conversion finished...")






#parser = GFFParser("/mnt/Data/GenomeData/Medicago/Annotation/Mt4.0v2_genes_20140818_1100.gff3")
#print("Features:", len(parser.features))

#parser = GFFParser("/mnt/Data/GenomeData/Nicotiana_benthamiana/annotation/Niben101_annotation.allfeatures.gff.gz")
#print("Features:", len(parser.features))