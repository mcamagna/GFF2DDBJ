
class FastaParser:
    
    def __init__(self, fasta_file_path):
        self.path = fasta_file_path
        self.parseFile()
        
        
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
        for i, line in enumerate(inp):
            if line[-1] == '\n':
                line = line[:-1]
            
            if line.startswith(">"):
                if i>0:
                    self.seqlens.append(current_seq_len)
                current_seq_len = 0
                
                line = line[1:].split(" ")[0]
                self.headers.append(line)
            else:
                current_seq_len+= len(line)
        self.seqlens.append(current_seq_len)
        inp.close()
        
        FastaParser.fasta_dict = dict()
        for h, length in zip(self.headers, self.seqlens):
            FastaParser.fasta_dict[h] = length
        
        
            
    def getFastaHeaders(self):
        return self.headers