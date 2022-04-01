# GFF2DDBJ
This (unofficial) tool will convert GFF3 files into DDBJ annotation files. Since DDBJ does not provide an official tool for such conversion, I created this tool to work with Braker2 generated GFF files, though it should work with most other GFF3 files as well.

## Disclaimer
Depending on how they were generated, GFF files differ substantially in their naming conventions, making an automated conversion into DDBJ annotations difficult. For this reason, this tool is not guaranteed to work correctly for all GFF3 files.
Also, the DDBJ annotation format is somewhat tedious and requires different mandatory parameters, such as special keywords, depending on the parameters provided in the header. The tool checks for various of these cases, but I have so far only implemented those scenarios of which I'm aware off.
  
**Use this tool at your own risk, and check the converted files thoroughly before submission**
DDBJ provides the UME tool, which will allow you to check whether your annotation files are correct. 

## Usage
Simply download this repository and execute the GFF2DDBJ.py file. There are no dependencies.

The most simple case is to simply run the tool using:
```
python GFF2DDBJ.py gff_file fasta_file 
```
<br>However, since DDBJ annotations begin with a mandatory set of meta-data, such as the submitters personal information, you will probably want to also provide your own data in the form of a text file. The format is identical to a DDBJ annotation file. See *example_header.txt* for further reference.
```
python GFF2DDBJ.py --header example_header.txt gff_file fasta_file 
```
<br>By default, only CDS features will be exported, but running this tool with the *--export_all* parameter will export all valid DDBJ features that can be found in the GFF file (mRNA, exons, intron, etc.).
```
python GFF2DDBJ.py --export_all gff_file fasta_file 
```
*Please note that exporting more features than absolutely necessary will make submission to DDBJ more difficult, since DDBJ enforces arbitrary rules, which are not (or not well) documented in their submission guidelines. While considerable efforts were made to satisfy their many criteria, this tool only incorporates rules that I was made aware off via email correspondence of my own WGS submission.*

<br><br>To see all available parameters, run
```
python GFF2DDBJ.py -h 
```

