#!/usr/bin/env python3
# Name: Justin Jang (jjang12)
# Group Members: Aster Lathbury(mlathbur), Kimberly Magpantay(klmagpan), Faiz Khan(faahkhan)
"""Analyze a genome"""

from sequenceAnalysis import NucParams, ProteinParam, FastAreader

class GenomeAnalyser:
    '''Display data from a Fasta file. Print sequence length, G-C content, and output statistics on 
    relative codon usage for each codon ordered by codon within each amino acid group. '''
    def __init__(self,fileName = 'testGenome.fa'):
        '''Take file and run it'''
        self.fastAreader = FastAreader()
        self.nucParams = NucParams()
        self.ProteinParam = ProteinParam()

        nucComp = NucParams.nucComposition()
        codonComp = NucParams.codonComposition()
        aaComp = NucParams.aaComposition()

    
    def analyser(self):
        for header, sequence in self.fastAreader.readFasta():
            NucParams.addSequence(sequence)


    
    def analyze(self):
        for header, sequence in NucParams.fastAReader.readFasta():
            NucParams.addSequence(sequence)

        '''get all the dictionaries we need for analysis'''
        nucComp = NucParams.nucComposition()
        codonComp = NucParams.codonComposition()
        aaComp = NucParams.aaComposition()

        '''extract data from dictionaries'''
        sequenceLength = NucParams.nucCount() / 1000000
        gcContent = 0
        if sequenceLength > 0:
            gcContent = ((nucComp['G'] + nucComp['C']) / sum(nucComp.values())) * 100

        '''output data'''
        print("sequence length = %0.2f Mb\n" % sequenceLength)
        print("GC content = %0.1f %% \n" % gcContent)

        for codon, aa in sorted(NucParams.rnaCodonTable.items()):
            codonFrequency = 0
            '''check if file has any aa in it'''
            if sequenceLength > 0:
                codonFrequency = (codonComp[codon] / aaComp[aa]) * 100
            '''print out information on codon bias'''
            print("%s : %s %5.1f (%6d)" % (codon, aa, codonFrequency, codonComp[codon]))

    def main():
        analyzer = GenomeAnalyser()
        analyzer.analyze()
    main()
    