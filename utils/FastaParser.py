import re

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
            import gzip
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
                if not line.startswith("\\\\"): #to ignore the completely unnecessary DDBJ fasta 'end flag'
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