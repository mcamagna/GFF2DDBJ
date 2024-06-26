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
from utils.features import TruncatedBothSidesFeature, CompoundFeature,\
    TruncatedFeature
from utils import GFFWriter

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
    #Added 2024/5 due to requirements from DDBJ
    parser.add_argument('--country', help="The country where the sample was collected.")
    parser.add_argument('--isolation_source', help="The isolation source of the sample.")
    parser.add_argument('--host', help="The host of the sample.")
    parser.add_argument('--collection_date', help="The collection date of the sample.")

    parser.add_argument('--mol_type', help="Type of molecule used in the sample. If not provided, you will be asked to choose the type if necessary.")
    parser.add_argument('--locus_tag_prefix', help="A prefix that is attached before each gene name. Must be 3-12 letters long and contain only alphanumeric characters. The first character should be a letter.")
    parser.add_argument('--export_all', action='store_true', help="Parses the GFF completely, but only writes the source and CDS features. For genome annotations this is typically sufficient and can avoid difficulties such as alternatative splicing, which is not handled well in DDBJ files.")
    parser.add_argument('--gene_as_note', action='store_true', help="By default, the gene name/id will be written as 'gene' qualifier into each feature belonging to that gene. Using this flag, each feature will instead be labeled with 'note gene ID' instead.")
    parser.add_argument('--intermediate_gff', help="Optional: Output path for the intermediate GFF file. During parsing of the GFF files, some changes to the information in the GFF file may need to be introduced to allow exporting the file. Writing this intermediate GFF file can be useful to track down sources of error.")
    
    #parser.print_help()
    
    
    args = parser.parse_args()
    INFILE = args.GFF
    FASTAFILE = args.FASTA
    OUTFILE = args.out
    Parameters.source_attributes['organism'] = args.organism
    Parameters.source_attributes['mol_type'] = args.mol_type
    Parameters.source_attributes['strain'] = args.strain
    Parameters.locus_attributes["locus_tag_prefix"] = args.locus_tag_prefix

    if args.country is not None:
        Parameters.source_attributes['country'] = args.country
    if args.collection_date is not None:
        Parameters.source_attributes['collection_date'] = args.collection_date
    if args.host is not None:
        Parameters.source_attributes['host'] = args.host
    if args.isolation_source is not None:
        Parameters.source_attributes["isolation_source"] = args.isolation_source

    Parameters.export_all = args.export_all
    Parameters.gene_as_note = args.gene_as_note
    Parameters.intermediate_gff = args.intermediate_gff
    
    
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
        #print("The COMMON header currently contains these values:")
        #Parameters.printCommonParameters()
    
    Parameters.askUserForRequiredParameters()
    
    ddbjwriter = DDBJWriter(OUTFILE)
    
    
    
    print("Parsing GFF file:", INFILE)
    gffparser = GFFParser(INFILE)
    print("Number of features found in GFF file:", len(gffparser.features))
    features = gffparser.features
    
    if Parameters.intermediate_gff is not None:
        GFFWriter.writeGFF(features)
    
    print("Parsing FASTA file")
    fastaParser = FastaParser(FASTAFILE)
    fasta_headers = fastaParser.getFastaHeaders()
    
    if len(fastaParser.assembly_gaps)>0:
        Parameters.askUserForAssemblyGapInfo()
        
    
    print("Converting features")
    fconverter = FeatureConverter()
    fconverter.convertFeatures(features)
    fconverter.addAssemblyGaps(features, fastaParser.assembly_gaps)
    

    features_to_translate = []
    for feature in features.values():
        if isinstance(feature, TruncatedBothSidesFeature) or (isinstance(feature, CompoundFeature) and isinstance(feature.members[0], TruncatedFeature) and isinstance(feature.members[-1], TruncatedFeature) and len(feature.members)>1):
            features_to_translate.append(feature)
    if len(features_to_translate)>0:
        print("Found coding sequences with missing start and stop codon. Guessing best reading frame... this may take a while.")        
        fastaParser.guessBestReadingFrame(features_to_translate)

    #Remove CDS entries that were flagged with an INVALID_CDS feature while guessing the best reading frame
    
    ddbjwriter.writeHeader()
    ddbjwriter.writeFeatures(features, fasta_headers)
    
    print("Conversion finished...")
    
    
if __name__ == "__main__":
    main()

