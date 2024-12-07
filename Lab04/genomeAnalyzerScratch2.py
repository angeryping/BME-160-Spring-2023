#!/usr/bin/env python3
# Name: Justin Jang (jjang12)
# Group Members: Aster Lathbury(mlathbur), Kimberly Magpantay(klmagpan), Faiz Khan(faahkhan)
from sequenceAnalysis import NucParams, FastAreader, ProteinParam

"""
This program prints statistics about the input file including
seq length, GC content, and codon usage percentages using sequenceAnalysis
This assignment taught me about imports, formatted printing, keyword self,
taking input using arguments, and dictionary loop operations.
Usage: python3 genomeAnalyzer.py SynechococcusGenome.fa 
"""


class genomeAnalyzer:
    """
    This program is the printer and last stage of math for
    argument fastAfile analysis.  It prints where the file is from,
    length of sequence, GC percentage, and codon usage
    """
    def __init__(self, filename='testGenome.fa'):
        """
        Takes in file from the directory and runs the analysis
        """

        self.analyser(filename)
    
    def gcContent(self, nucParams):
        """
        This function returns gc content percentage
        """
        nucComp = nucParams.nucComposition()
        GC = nucComp['G'] + nucComp['C']
        GC = GC / nucParams.nucCount() * 100
        return GC

    def sequenceLength (self, myRead):
        """
        This program takes in a fasta file, runs it through NucParams
        and returns length divided by a million to return megabases
        """
        nucParams = NucParams('')
        for head, seq in myRead.readFasta():
            nucParams.addSequence(seq)

        length = nucParams.nucCount()/ 1000000

        return length, nucParams

    def analyser(self, fastAfile):
        """
        This method does all the printing and calculates codon usage
        """
        # print("Reading from {}...".format(fastAfile))
        self.myReader = FastAreader(fastAfile)

        len, nucParams = self.sequenceLength(self.myReader)
        print((f"sequence length = {len:.2f} Mb\n"))

        self.GC = self.gcContent(nucParams)
        print(f"GC content = {self.GC:.1f} %\n")

        codonCount = nucParams.codonComposition()
        totalAAs = 0

        # Alphabatizes and loops through codonCount keys
        for codons in sorted(codonCount):
            c_count = codonCount[codons]
            totalAAs += c_count
            one_letter = nucParams.rnaCodonTable[codons]
            aminoCount = nucParams.aaComp[one_letter]
            #calculates percentage
            if c_count != 0:
                finalPercentage = (c_count/aminoCount) *100
            else:
                finalPercentage = (c_count/1) *100
            print('{} : {} {:5.1f}% ({:6d})'.format(codons, one_letter, finalPercentage, c_count))
        #print("Total number of amino acids counted is "+str(totalAAs))


    
genomeAnalyzer()