#!/usr/bin/env python3
# Name: Justin Jang (jjang12)
# Group Members: Aster Lathbury(mlathbur), Kimberly Magpantay(klmagpan), Faiz Khan(faahkhan)
"""Analyze a genome from a FastaAfile."""

from sequenceAnalysis import NucParams, ProteinParam, FastAreader
import sys

class GenomeAnalyser:
    '''Display data from a Fasta file, in this case "testGenome.fa". Print sequence length, G-C content, and output statistics on 
    relative codon usage for each codon ordered by codon within each amino acid group. '''

    def __init__(self,fileName = 'testGenome.fa'):
        '''Constructor for GenomeAnalyser class.Take file and run the analyser method on it.'''
        #useful dictionaries might use?
        #add more stuff if needed
        self.fastAreader = FastAreader()  
        self.nucParams = NucParams()

        if fileName == None:  #in case fileName is none
            fileName = sys.stdin
        else:
            self.analyser(fileName)
    
    
    def sequenceLength(self,read):
        '''Run NucParams on the fastA file and return megabases'''
        nucParams = NucParams()  # easier way to call NucParams methods
        for head,seq in read.readFasta():  # note to self: takes readFasta from FastAreader in sequenceAnalysis
            nucParams.addSequence(seq)

        length = nucParams.nucCount() / 1000000  # return megabases by dividing by 1,000,000
        return length, nucParams
    
    def gcContent(self, nucParams):
        '''Calculate the percentage of G and C nucleotides in sequence from the file.'''
        nucComp = nucParams.nucComposition() #see NucParams in sequenceAnalysis
        gc = nucComp['G'] + nucComp['C']
        gc = gc / nucParams.nucCount() * 100  #see NucParams in sequenceAnalysis

        return gc

    def analyser(self, fastAfile):
        self.read = FastAreader(fastAfile)

       
        len, nucParams = self.sequenceLength(self.read)  # call sequenceLength with self.read as the argument, return len = length of sequence, nucParams = instance of NucParams with compositions stats
        print(f"sequence length = {len:.2f} Mb\n")

        
        self.GC = self.gcContent(nucParams)  # call gcContent with nucParams as the argument, which assigns GC content to the variable self.GC
        print(f"GC content = {self.GC:.1f} %\n")

        codonCount = nucParams.codonComposition() # store codon composition stats in codonCount
        

        for codon in sorted(codonCount.keys(), key=lambda c: (nucParams.rnaCodonTable[c], c)):  # csort alphabetically first by 1-letter amino acid code, THEN by 3-letter amino acid code.
        #is this a troll way of commenting idk                                                                                        # c represents each codon in codonCount
            cCount = codonCount[codon] # count of current codon, store in cCount
            singLetter = nucParams.rnaCodonTable[codon] # store single-letter AA in singLetter
            aminoCount = nucParams.aaComp[singLetter] # store count of AA in aminoCount

            # calculate frequencies
            if cCount > 0:
                finalFreq = (cCount/aminoCount) *100
                print(f'{codon} : {singLetter}  {finalFreq:5.1f}% ({cCount:6d})')
            else:
                finalFreq = (cCount) *100
                print(f'{codon} : {singLetter}  {finalFreq:5.1f}% ({cCount:6d})')

        #print(codonCount)  run this to test raw count before showing frequencies
        #test for redundancy later
    
GenomeAnalyser()

#compare results with test.out    