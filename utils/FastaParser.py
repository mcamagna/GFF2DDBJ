import re, gzip
from utils.features import CompoundFeature

class FastaParser:
    
    def __init__(self, fasta_file_path):
        self.path = fasta_file_path
        self.gap_regex = re.compile("(N|n)+")
        self.assembly_gaps = dict()
        self.parseFile()
    
    def findGaps(self, fasta_entry_name, seq):    
        gap_list = None
        for gap in self.gap_regex.finditer(seq):
            if gap_list is None: #this is to assure that we only fetch the list once for performance reasons
                gap_list = self.assembly_gaps.get(fasta_entry_name)
                if gap_list is None:
                    gap_list = []
                    self.assembly_gaps[fasta_entry_name] = gap_list
            gap_list.append((gap.start(), gap.end()))
    
        
    def parseFile(self):
        inp = None
        if self.path.lower().endswith("gz") or self.path.lower().endswith("gzip"):
            inp = gzip.open(self.path, 'rt')
        else:
            inp = open(self.path, 'rt')
            
        self.headers = []
        self.seqlens = []
        
        current_seq_len = 0
        current_lines = []
        
        for i, line in enumerate(inp):
            if line[-1] == '\n':
                line = line[:-1]
            
            if line.startswith(">"):
                line = line[1:].split(" ")[0]
                
                if i>0:
                    current_seq = "".join(current_lines)
                    current_lines.clear()
                    self.findGaps(self.headers[-1], current_seq)
                    current_seq = None
                    self.seqlens.append(current_seq_len)
                current_seq_len = 0
                self.headers.append(line)
                
            else:
                line = line.replace('\\', '').replace('//','') #to ignore the completely unnecessary DDBJ fasta 'end flag'
                current_seq_len+= len(line)
                current_lines.append(line)
            
        #Last fasta entry needs to be added        
        current_seq = "".join(current_lines)
        current_lines.clear()
        self.findGaps(self.headers[-1], current_seq)
        current_seq = None
        self.seqlens.append(current_seq_len)
        inp.close()
        
        FastaParser.fasta_dict = dict()
        for h, length in zip(self.headers, self.seqlens):
            FastaParser.fasta_dict[h] = length
        
    def getFastaHeaders(self):
        return self.headers
    
    
    
    
    def reverseComplement(self, sequence):
        sequence = sequence.upper()
        rc = ""
        for s in reversed(sequence):
            if s=="A":
                rc+="T"
            elif s == "T":
                rc+= "A"
            elif s == "C":
                rc+= "G"
            elif s == "G":
                rc+= "C"
        return rc
    
    def extractSequence(self, feature, genomeseq):
        extracted = ""
        if isinstance(feature, CompoundFeature):
            for m in feature.members:
                extracted += genomeseq[m.start-1: m.end] #python starts position at 0
        else:
            extracted += genomeseq[feature.start-1 : feature.end]
        
        extracted = extracted.upper()
       
        if feature.strand == "-":
            extracted = self.reverseComplement(extracted)
            
        return extracted
    
    def evaluateReadingFrames(self, cds_seq):
        
        pattern = re.compile("TGA|TAA|TAG")
        best_frame = -1
        best_stopcodons = 999
        
        for frame in [0,1,2]:
            current_stopcodons = 0
            orf = cds_seq[frame:]
            x = 0
            while x<len(orf)-3:
                codon = orf[x:x+3]
                if len(pattern.findall(codon))>0:
                    current_stopcodons+=1
                x+=3
            
            if current_stopcodons<= best_stopcodons:
                best_stopcodons = current_stopcodons
                best_frame = frame
                
        return best_frame
        
        
    def translate_sequence(self, dna_seq):
        codontab = {
                    'TCA': 'S','TCC': 'S','TCG': 'S','TCT': 'S','TTC': 'F','TTT': 'F','TTA': 'L','TTG': 'L','TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C',
                    'TGT': 'C','TGA': '*','TGG': 'W','CTA': 'L','CTC': 'L','CTG': 'L','CTT': 'L','CCA': 'P','CCC': 'P','CCG': 'P','CCT': 'P','CAC': 'H','CAT': 'H','CAA': 'Q','CAG': 'Q','CGA': 'R',
                    'CGC': 'R','CGG': 'R','CGT': 'R','ATA': 'I','ATC': 'I','ATT': 'I','ATG': 'M','ACA': 'T','ACC': 'T','ACG': 'T','ACT': 'T','AAC': 'N','AAT': 'N','AAA': 'K','AAG': 'K','AGC': 'S',
                    'AGT': 'S','AGA': 'R','AGG': 'R','GTA': 'V','GTC': 'V','GTG': 'V','GTT': 'V','GCA': 'A','GCC': 'A','GCG': 'A','GCT': 'A','GAC': 'D','GAT': 'D','GAA': 'E','GAG': 'E','GGA': 'G',
                    'GGC': 'G','GGG': 'G','GGT': 'G' }
                    
        prot_seq = ""
        for i in range(0, len(dna_seq)-2, 3):
            codon = dna_seq[i]+ dna_seq[i+1]+dna_seq[i+2]
            aa = codontab.get(codon.upper())
            if aa is None:
                print(f"WARNING: Could not translate codon {codon}.")
            else:
                prot_seq += aa
        print(prot_seq)
        return prot_seq


    def guessBestReadingFrame(self, ddbj_features):
        """For features where both start and end positions are unkown, we need to obtain
        the DNA sequence, and check all three codon offsets for whether we get stop codons within
        the sequence"""
        
        fasta_headers_of_interest = set()
        for feature in ddbj_features:
            fasta_headers_of_interest.add(feature.seqid)
        
        inp = None
        if self.path.lower().endswith("gz") or self.path.lower().endswith("gzip"):
            
            inp = gzip.open(self.path, 'rt')
        else:
            inp = open(self.path, 'rt')
        
        currentHeader = ""
        currentSeq = ""
        foundSequenceOfInterest = False
        for line in inp:
            if line.startswith("\\\\") or line.startswith('//'):
                continue
            
            if line[-1] == '\n':
                line = line[:-1]  
                
            
            if line.startswith(">"):
                if foundSequenceOfInterest:#if the current sequence needs to be parsed
                    for i, feature in enumerate(ddbj_features):
                        if feature.seqid == currentHeader:
                            extracted = self.extractSequence(feature, currentSeq)
                            best_frame = self.evaluateReadingFrames(extracted)
                            feature.attributes["codon_start"] = str(best_frame+1)
                            #TODO: Check if contains stop codon and adjust CDS range otherwise
                            print(feature.attributes)
                            #print(extracted)
                            prot_seq = self.translate_sequence(extracted) 
                            print(prot_seq)
                            #if "*" in prot_seq[:-1]:
                            if "*" in prot_seq:
                                print("WARNING: Found stop codon within translated CDS sequence.")
                                feature.attributes['INVALID_CDS'] = "INVALID_CDS"
                            
                currentHeader = line[1:].split(" ")[0]
                currentSeq=""
                foundSequenceOfInterest = currentHeader in fasta_headers_of_interest
                
            elif not foundSequenceOfInterest:
                continue 
            else:
                currentSeq+= line
                
        inp.close()
        