
class FastaParser:
    
    def __init__(self, fasta_file_path):
        self.path = fasta_file_path
        
        
    def getFastaHeaders(self):
        inp = None
        if self.path.lower().endswith("gz") or self.path.lower().endswith("gzip"):
            import gzip
            inp = gzip.open(self.path, 'rt')
        else:
            inp = open(self.path, 'rt')
            
        headers = []
        for line in inp:
            if line.startswith(">"):
                line = line[1:].replace("\n","").split(" ")[0]
                headers.append(line)
        inp.close()
        return headers