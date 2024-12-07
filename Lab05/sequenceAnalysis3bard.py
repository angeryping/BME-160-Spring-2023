import sys
class FastAreader :
   
    def __init__ (self, fname=None):
        
        self.fname = fname
            
    def doOpen (self):
      
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
      
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
          
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



class OrfFinder():  # place in sequenceAnalysis when complete
    """Turn FASTA sequence into lists of lists of Open Reading Frames (ORFs)"""
    stopCodons = ['TGA', 'TAG', 'TAA']   # list of stop codons
    startCodons = ['ATG']   #list of start codons
    nucComplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}  # dict of valid nucleotides

    def __init__(self, genome):
        """Set up a list of lists of ORFs """
        self.genome = genome
        self.orfs = []  # store ORFs frame, start, stop, and length

    def keepOrf(self, frame, start, stop, length):
        """ Saves ORF info
        """
        self.orfs.append([frame, start, stop, length])

    def findOrfs(self):
        """Find ORFs on a 3'-5'("top") strand. Return it as a list of ORFs."""
        startPositions= []  #store start positions
        foundStart = False  # if a codon is found, will set to True for use in conditionals
        foundCodon = False 

        for frame in range(0,3):
            starts = []
            foundStart = False
            foundCodon = False
            for p in range(frame,len(self.genome),3):
                codon = self.genome[p: p + 3]
                if codon == 'ATG': #in OrfFinder.startCodons:
                    startPositions.append(p)
                    foundStart = True
                    foundCodon = True
                
                if codon in self.stopCodons and foundStart:
                    stop = p + 3
                    start = p
                    self.keepOrf((frame % 3) + 1, startPositions[0] + 1 - frame, p + 3, stop - start + 1)
                    startPositions = []
                    foundStart = False
                    foundCodon = True
                
            if foundStart:
                stop = p + 3
                start = p
                self.keepOrf((frame % 3) + 1, startPositions[0] + 1, len(self.genome), stop - start + 1)

        return self.orfs
    
    def reverseComplement(self):
        """ Create the reversed and complimentary strand of DNA
        """
        return ''.join([self.nucComplement[base] for base in self.genome[::-1]])

    def revOrfs(self):
        """Create the reverse strand 5'-3' strand. """
        reverseComp = self.reverseComplement()
        start_positions = []
        foundStart = 0
        foundCodon = 0

        for frame in range(0, 3):
            foundStart = 0 
            foundCodon = 0
            start_positions = []
            for p in range(frame, len(reverseComp), 3):
                codon = reverseComp[p:p + 3]
                stop = p + 3
                start = p

                if codon == 'ATG':
                    start_positions.append(p)
                    foundStart = 1
                    foundCodon = 1
                if codon in self.stopCodons and foundStart:
                    self.keepOrf(-1 * ((frame % 3) + 1), len(reverseComp) - start_positions[0], len(reverseComp) - (p+2), stop - start + 1)
                    start_positions = []
                    foundStart = False
                    foundCodon = True
                if not foundCodon and codon in self.stopCodons:
                    self.keepOrf(-1 * ((frame % 3) + 1), len(reverseComp) - p - 2, len(reverseComp), stop - start + 1)
                    start_positions = []
                    foundCodon = True
            if foundStart:
                self.keepOrf(-1 * ((frame % 3) + 1), start_positions[0] + 1, 1, stop - start + 1)
            
        return self.orfs