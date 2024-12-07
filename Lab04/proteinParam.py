#!/usr/bin/env python3
# Name: Justin Jang (jjang12)
# Group Members: Aster Lathbury(mlathbur), Kimberly Magpantay(klmagpan), Faiz Khan(faahkhan)

'''
Take a protein string and calculate the physical-chemical properties of a protein sequence. Return:
- the number of amino acids and total molecular weight,
- Molar extinction coefficient and Mass extinction coefficient,
- theoretical isoelectric point (pI)
- amino acid composition

Input: 
VLSPADKTNVKAAW

Output: 
Number of Amino Acids: 14
Molecular Weight: 1499.7
molar Extinction coefficient: 5500.00
mass Extinction coefficient: 3.67
Theoretical pI: 9.88
Amino acid composition:
A = 21.43%
C = 0.00%
D = 7.14%
E = 0.00%
F = 0.00%
G = 0.00%
H = 0.00%
I = 0.00%
K = 14.29%
L = 7.14%
M = 0.00%
N = 7.14%
P = 7.14%
Q = 0.00%
R = 0.00%
S = 7.14%
T = 7.14%
V = 14.29%
W = 7.14%
Y = 0.00%
'''

class ProteinParam :
    """Define a class named ProteinParam. Defines dictionaries and variables in class attributes."""
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    #dictionary of amino acids and corresponding molecular weights
    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    #molecular weight water
    mwH2O = 18.015
    #AAs and their absorbance
    aa2abs280 = {'Y':1490, 'W': 5500, 'C': 125}

    #charged AAs and their pKa vals
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    #pKas of N- and C- terminus
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, proteinStr: str):
        """Initialize method for taking protein sequence as an inputted string"""

        #empty list to hold individual AAs in protein sequence
        splitAminos = []

        #list allowedAminos is the AAs from aa2mw
        allowedAminos = self.aa2mw.keys()

        for char in proteinStr.upper():  # input to upper
            #iterate thru each char in input str, 
            #if char in input, then append to splitAminos list.
            if char in allowedAminos:
                splitAminos.append(char)

        #join strings in splitAminos, then separate them into individual AAs in the list.
        list = ''.join(splitAminos)
        # joins strings with no separator adding: ''
        # splits joined 'protein' into a char array, assigns to l

        ## Join the list of amino acids back into a single string, and convert to uppercase.
        self.newProteinStr = ''.join(list).upper()
        # initialize newProteinStr to joined l and makes uppercase
        
        # Create an empty dictionary to hold the amino acid composition of the protein.
        self.aaComp = {}
        
        # creates an empty dictionary object to represent # of each amino Acid
        for aminoAcid in self.aa2mw.keys():
            # initializes aminoAcid, loops through aa2mw.keys
            self.aaComp[aminoAcid] = self.newProteinStr.count(aminoAcid)
            # assigns aminoAcid to the key, counts self.newProteinStr for the amino acid letter

    def aaCount (self):
        """Return a single integer count of valid amino acid characters found."""
        return (len(self.newProteinStr))
        # returns len of input self's newProteinStr object

    def _charge_ (self,pH):
        """ Calculate the net charge of a protein based on its amino acid composition and the pH of the surrounding environment. Use for pI method."""   
        #ask TA questions about _charge_
        posCharge = sum(self.aaComp.get(aa, 0) * pow(10, ProteinParam.aa2chargePos.get(aa, 0)) / (pow(10, ProteinParam.aa2chargePos.get(aa, 0)) + pow(10, pH)) for aa in ['R', 'K', 'H'])  # calculate positive charge from equation
        negCharge = sum(self.aaComp.get(aa, 0) * pow(10, pH) / (pow(10, ProteinParam.aa2chargeNeg.get(aa, 0)) + pow(10, pH)) for aa in ['D', 'E', 'C', 'Y'])  # calculate negative charge from equation
        
        nTermChrg = pow(10, ProteinParam.aaNterm) / (pow(10,ProteinParam.aaNterm) + pow(10,pH))
        cTermChrg = pow(10,pH) / (pow(10,ProteinParam.aaCterm) + pow(10,pH))

        netCharge = posCharge - negCharge + nTermChrg - cTermChrg
        return netCharge

    def pI(self):
        """Calculate the  particular pH that yields a neutral net charge"""
        bestCharge = 100000000
        bestpH = 0
        for pH100 in range (0, 1400+1): # want this down to two precision points, so much put it from 1 to 1400 not 1 to 14
            pH = (pH100)/100 #gets the two decimal numbers (gives actual pH)
            currentCharge = self._charge_(pH) #runs pH through the charge method
            keepCharge = abs(currentCharge) #call charge method, gets absolute value 
            
            if keepCharge < bestCharge: # wants the closest value to 0, so keeps it if it is lower than previous
                bestCharge = keepCharge #defines keepCharge as new lower value, runs it through the loop again 
                bestpH =  pH
        return bestpH   
        
    

    def aaComposition (self) :
        """return a dictionary keyed by single letter Amino acid code, and having associated values that are the counts of those amino acids in the sequence."""
        return self.aaComp #returns the dictionary setup at the beginning

    def molarExtinction (self):
        """Return how much light a protein absorbs at a certain wavelength"""
        return (self.newProteinStr.count('Y') * self.aa2abs280['Y']
                     + self.newProteinStr.count('W') * self.aa2abs280['W']
                     + self.newProteinStr.count('C') * self.aa2abs280['C'])
        

    def massExtinction (self):
        """
        Calculate light absorbance divided by molecular weight of the input protein sequence.
        """
        myMW =  self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0
    

    def molecularWeight (self):
        """Return molecular weight of the protein sequence."""
        waterWeight = self.mwH2O * (self.aaCount() - 1)
        aaWeight = sum(count * self.aa2mw[aa] for aa, count in self.aaComp.items())
        molecularWeight = aaWeight - waterWeight
        return molecularWeight
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
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present 
        
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?')

if __name__ == "__main__":
    main()