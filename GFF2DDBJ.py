'''
@author: Maurizio Camagna
'''
import sys, os
from utils.GFFParser import GFFParser
from utils.DDBJWriter import DDBJWriter
from utils.FeatureConverter import FeatureConverter
from utils.FastaParser import FastaParser
from utils.Parameters import Parameters
import argparse

def checkFilepaths(filepaths):
    for path in filepaths:
        
        if not os.path.exists(path):
            print("ERROR: No file found at", path)
            sys.exit(1)


def main():
    Parameters.init()
    parser = argparse.ArgumentParser(description='A tool to help you convert GFF3 files into DDBJ annotation files.')
    parser.add_argument('GFF', help='Path to a GFF3 file (can be gzipped).')
    parser.add_argument('FASTA', help='Path to a FASTA file (can be gzipped).')
    parser.add_argument('--out', help="Optional: Location where the DDBJ annotation will be stored. If nothing is provided, the annotation file will be stored in the same location as the GFF3 file.")
    parser.add_argument('--header', help="Optional: Location of the text file specifying the values for the DDBJ header. Check example_header.txt for more information.")
    parser.add_argument('--organism', help="Optional: Scientific name of the organism.")
    parser.add_argument('--strain', help="Name of the strain.")
    parser.add_argument('--mol_type', help="Type of molecule used in the sample. If not provided, you will be asked to choose the type if necessary.")
    parser.add_argument('--locus_tag_prefix', help="A prefix that is attached before each gene name. Must be 3-12 letters long and contain only alphanumeric characters. The first character should be a letter.")
    parser.add_argument('--no_sorting', action='store_true', help="By default, the features of each source entry will be sorted by position. Use this flag to disable this behaviour.")
    
    #parser.print_help()
    
    
    args = parser.parse_args()
    INFILE = args.GFF
    FASTAFILE = args.FASTA
    OUTFILE = args.out
    Parameters.source_attributes['organism'] = args.organism
    Parameters.source_attributes['mol_type'] = args.mol_type
    Parameters.source_attributes['strain'] = args.strain
    Parameters.locus_attributes["locus_tag_prefix"] = args.locus_tag_prefix
    Parameters.sort_features = not args.no_sorting
    
    if OUTFILE is None:
        OUTFILE = INFILE.replace(".gff3", "").replace(".GFF3", "").replace(".gff", "").replace(".GFF", '')
        OUTFILE += ".ann"
        print("Annotation will be written to:", OUTFILE)
    
    HEADERFILE = args.header
    if HEADERFILE is None:
        print("Warning: No header file was provided. Make sure to manually add the header after the conversion.")
        checkFilepaths([INFILE, FASTAFILE])
    else:
        checkFilepaths([INFILE, FASTAFILE, HEADERFILE])
        Parameters.parseHeaderFile(HEADERFILE)
        print("The COMMON header currently contains these values:")
        Parameters.printCommonParameters()
    Parameters.askUserForRequiredParameters()
    
    ddbjwriter = DDBJWriter(OUTFILE)
    
    
    
    print("Parsing GFF file:", INFILE)
    gffparser = GFFParser(INFILE)
    print("Number of features found in GFF file:", len(gffparser.features))
    features = gffparser.features
    
    print("Parsing FASTA file")
    fastaParser = FastaParser(FASTAFILE)
    fasta_headers = fastaParser.getFastaHeaders()
    
    if len(fastaParser.assembly_gaps)>0:
        Parameters.askUserForAssemblyGapInfo()
        
    
    print("Converting features")
    fconverter = FeatureConverter()
    features = fconverter.convertFeatures(features)
    features = fconverter.addAssemblyGaps(features, fastaParser.assembly_gaps)
    
    ddbjwriter.writeHeader()
    ddbjwriter.writeFeatures(features, fasta_headers)
    
    print("Conversion finished...")
    
    
if __name__ == "__main__":
    main()

