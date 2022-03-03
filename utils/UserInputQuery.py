'''
This will ask the user for additional attributes and suggest correct values
@author: Maurizio Camagna
'''

import sys

class UserInputQuery:
    def __init__(self):
        self.attributes = dict()
    
        
        
    
    def askUserForMolType(self):
        print("\nI couldn't find the mandatory value 'mol_type' in the header.")
        print("Please select a mol_type from the list:")
        mol_types = {"1":"genomic DNA", 
                      "2":"genomic RNA", 
                      "3":"mRNA", 
                      "4":"tRNA", 
                      "5":"rRNA", 
                      "6":"other RNA", 
                      "7":"other DNA", 
                      "8":"transcribed RNA", 
                      "9":"viral cRNA", 
                      "10":"unassigned DNA", 
                      "11":"unassigned RNA",
                      "0": 'EXIT'}
        for key in mol_types.keys():
            print("  ["+str(key)+"] "+mol_types[key])
        
        inp = input("Select mol_type number: ")
        selection = mol_types.get(inp)
        if selection is not None and selection !='EXIT':
            return selection
        else:
            print("ERROR: Invalid selection... exiting")
            sys.exit(128)
        
    def askUserForOrganism(self):
        print("\nI couldn't find the mandatory value 'organism' in the header.")
        inp = input("Please enter the name of the organism: ")
        return inp

        
  