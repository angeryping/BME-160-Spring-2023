class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        upperProt = str.upper(protein)
        splitAminos = []
        allowedAminos = self.aa2mw.keys()
        for char in upperProt:
            if char in allowedAminos:
                splitAminos.append(char)

        list = ''.join(splitAminos).split()
        # joins strings with no separator adding: ''
        # splits joined 'protein' into a char array, assigns to l
        self.newProteinStr = ''.join(list).upper()
        # initialize newProteinStr to joined l and makes uppercase

        self.aaComp = {}
        # creates an empty dictionary object to represent # of each amino Acid
        for aminoAcid in self.aa2mw.keys():
            # initializes aminoAcid, loops through aa2mw.keys
            self.aaComp[aminoAcid] = self.newProteinStr.count(aminoAcid)
            # assigns aminoAcid to the key, counts self.newProteinStr for the amino acid letter

    def aaCount (self):
        return (len(self.newProteinStr))
        # returns len of input self's newProteinStr object

    def _charge_ (self,pH):
        posCharge = 0
        negCharge = 0.0 
        
        
        #for aa, pka in self.aa2chargePos.items():
        posCharge = sum(self.aaComposition.get(aa, 0) * pow(10, ProteinParam.aa2chargePos.get(aa, 0)) / (pow(10, ProteinParam.aa2chargePos.get(aa, 0)) + pow(10, pH)) for aa in ['R', 'K', 'H', ProteinParam.aaNterm])
        #for aa, pka in self.aa2chargeNeg.items():
        negCharge = sum(self.aaComposition.get(aa, 0) * pow(10, pH) / (pow(10, ProteinParam.aa2chargeNeg.get(aa, 0)) + pow(10, pH)) for aa in ['D', 'E', 'C', 'Y', ProteinParam.aaCterm])
      

    def pI(self):
      bestCharge = 100000000
      for pH100 in range(0,1400+1):
        pH = pH100 / 100 # pH is equal to 14
        thisCharge = abs(charge(pH))
        if thisCharge < bestCharge :
          bestCharge = thisCharge
      return bestCharge

    def aaComposition (self) :
        return self.aaComp #returns the dictionary setup at the beginning

    def molarExtinction (self):
        return (self.newProteinStr.count('Y') * self.aa2abs280['Y']
                     + self.newProteinStr.count('W') * self.aa2abs280['W']
                     + self.newProteinStr.count('C') * self.aa2abs280['C'])

    def massExtinction (self):
        """
        massExtiction does not cause the die off of many species as the name suggests.
        It calculates light absorbance divided by molecular weight of the input protein
        """
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        return self.molecularWeight == sum([ProteinParam.aa2mw[aa] for aa in self.newProteinStr]) + ProteinParam.mwH2O * (len(self.newProteinStr) - 1)

# Please do not modify any of the following.  This will produce a standard output that can be parsed
    
import sys
def main():
    inString = input('protein sequence?: ')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        #print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()
      
    #VLSPADKTNVKAAW
