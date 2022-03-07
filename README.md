# GFF2DDBJ
This (unofficial) tool will convert GFF3 files into DDBJ annotation files. Since DDBJ does not provide an official tool for such conversion, I created this tool to work with Braker2 generated GFF files, though it should work with most other GFF3 files as well.

## Disclaimer
Depending on how they were generated, GFF files differ substantially in their naming conventions, making an automated conversion into DDBJ annotations difficult. For this reason, this tool is not guaranteed to work correctly for all GFF3 files.
Also, the DDBJ annotation format is somewhat tedious and requires different mandatory parameters, such as special keywords, depending on the parameters provided in the header. The tool checks for various of these cases, but I have so far only implemented those scenarios of which I'm aware of.
  
**Use this tool at your own risk, and check the converted files thoroughly before submission**
DDBJ provides the UME tool, which will allow you to check whether your annotation files are correct. 

## Usage
Simply download this repository and execute the GFF2DDBJ.py file. There are no dependencies.

The most simple case is to simply run the tool using:

python GFF2DDBJ.py *path_to_gff_file* *path_to_fasta_file* 
