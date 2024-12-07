#!/usr/bin/env python3
# Name: Justin Jang (jjang12)
# Group Members: Aster Lathbury(mlathbur), Kimberly Magpantay(klmagpan), Faiz Khan(faahkhan)
"""
Analyze a genome. R


Methods: 
NucParams: Return counts of codons and their translated amino acids.
ProteinParam: Return aa count, pI, molar/mass extinction, and molecular weight
FastAreader: Read FastA files so we can take objects from the file to analyze in genomeAnalyzer




"""




import sys
class NucParams:
    """Track counts of codons and their translated amino acids.
     nucleotide > codon
     Methods:
     _init_(self)
     addSequence(self,sequence)
     aaComposition(self)  (from NB)
     nucComposition(self)  (from NB)
     codonComposition(self)  (from NB)
     nucCount(self)  (from NB)
       """
    rnaCodonTable = {
    # RNA codon table
    # U  
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        '''Set up dictionaries, set value to zero. Add a tally if sequence contains an amino acid.'''
        
        
        #make dictionaries with keys of sequence composition, values counts set to 0
        self.nucComp = {'G': 0,'U': 0,'A': 0,'C': 0,'T': 0,'N': 0}  # nucleotide composition
        self.aaComp = {aa: 0 for aa in NucParams.rnaCodonTable.values()}  # amino acid composition, key is amino acid and value is set to 0
        self.codonComp = {codon: 0 for codon in NucParams.rnaCodonTable.keys()} #codon sequence composition, key is codon and value is set to 0
        self.addSequence(inString)  # call addSequence method with the inputted sequence
               
        
    def addSequence (self, sequence):
        '''Clean the input sequence and increment how many of each RNA codon and its amino acid symbols in a list.
         Store in dictionaries self.codonComp and self.aaComp'''
        
        sequence = ''.join(sequence.split()).upper()  # remove whitespace, convert to uppercase

        for base in self.nucComp:  # update nucleotide composition (nucComp) with # of each nucleotide in sequence
            self.nucComp[base] += sequence.count(base)  #
        
       
        cleanedCodons = []  # create empty list to hold 3-base codons
        for i in range(0, len(sequence), 3):  # iterate over sequence inputted string in increments of 3 nucleotides
            cleanedCodons.append(sequence[i:i+3])  # place at position 0, 3, 6, 9 ...
        

        for codonSeq in cleanedCodons:  # iterate over cleanedCodons[] list
            codonSeq = codonSeq.replace('T','U')  # change from DNA to RNA format b/c we are matching sequence with codonComp which is RNA format
            if 'N'in codonSeq:  # skip N
                continue
            if codonSeq in self.codonComp:  # compare sequence with codon composition dictionary and searches for a match
                self.codonComp[codonSeq] += 1  # increment codon sequence value by 1 IF there is a match
                
                aa = self.rnaCodonTable[codonSeq] # 
                self.aaComp[aa] += 1 #increment amino acid value by 1
        pass
    def aaComposition(self):
        '''Return a composition of amino acids'''
        return self.aaComp
    def nucComposition(self):
        '''Return a compostion of nucleotides'''
        return self.nucComp
    def codonComposition(self):
        '''Return a compostion of 3-letter codons'''
        return self.codonComp
    def nucCount(self):
        '''Return a count of each kind of nucleotide'''
        return sum(self.nucComp.values())

class ProteinParam :
    """Calculate statistics of a given protein sequence. Take a protein string and calculate the physical-chemical properties of a protein sequence. 
Return:
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
Y = 0.00%"""
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


import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

