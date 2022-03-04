'''
This will ask the user for additional attributes and suggest correct values
@author: Maurizio Camagna
'''

import sys

class UserInputQuery:
    def __init__(self):
        self.attributes = dict()
    
        
        
    def askForWGS(self):    
        print("\nThe DATATYPE states WGS. One of the following selections must be chosen:")
 
        responses = {"1":"STANDARD_DRAFT", 
                      "2":"HIGH_QUALITY_DRAFT", 
                      "3":"IMPROVED_HIGH_QUALITY_DRAFT", 
                      "4":"NON_CONTIGUOUS_FINISHED",
                      "0": 'EXIT'}
        for key in responses.keys():
            print("  ["+str(key)+"] "+responses[key])
        inp = input("Select number: ")
        selection = responses.get(inp)
        if selection is not None and selection !='EXIT':
            return selection
        else:
            print("ERROR: Invalid selection... exiting")
            sys.exit(128)
            
            
    def askForTPA(self):    
        print("\nThe DATATYPE states TPA. One of the following selections must be chosen:")
 
        responses = {"1":"TPA:inferential", 
                      "2":"TPA:experimental", 
                      "0": 'EXIT'}
        for key in responses.keys():
            print("  ["+str(key)+"] "+responses[key])
        inp = input("Select number: ")
        selection = responses.get(inp)
        if selection is not None and selection !='EXIT':
            return selection
        else:
            print("ERROR: Invalid selection... exiting")
            sys.exit(128)
            
            
    def askForEST(self):    
        print("\nYou specified this as EST submission. One of the following selections must be chosen:")
 
        responses = {"1":"5’-end sequence (5’-EST)", 
                      "2":"3’-end sequence (3’-EST)", 
                      "0": 'EXIT'}
        for key in responses.keys():
            print("  ["+str(key)+"] "+responses[key])
        inp = input("Select number: ")
        selection = responses.get(inp)
        if selection is not None and selection !='EXIT':
            return selection
        else:
            print("ERROR: Invalid selection... exiting")
            sys.exit(128)
    
    
    def askUserForMolType(self):
        print("\nI couldn't find the mandatory value 'mol_type' in the header.")
        print("Please select a mol_type from the list:")
        responses = {"1":"genomic DNA", 
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
        for key in responses.keys():
            print("  ["+str(key)+"] "+responses[key])
        
        inp = input("Select mol_type number: ")
        selection = responses.get(inp)
        if selection is not None and selection !='EXIT':
            return selection
        else:
            print("ERROR: Invalid selection... exiting")
            sys.exit(128)
        
    def askUserForOrganism(self):
        print("\nI couldn't find the mandatory value 'organism' in the header.")
        inp = input("Please enter the name of the organism: ")
        return inp


    def askUserForStrain(self):
        print("")
        inp = input("Please enter the name of the strain: ")
        return inp
  