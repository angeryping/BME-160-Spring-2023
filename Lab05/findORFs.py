#!/usr/bin/env python3
# Name: Justin Jang (jjang12)
# Group Members: Aster Lathbury(mlathbur), Kimberly Magpantay(klmagpan), Faiz Khan(faahkhan)
"""
findORFs
Analyze a FASTA-formatted file containing a sequence of DNA and find the ORFs (start and stop codons).
Design notes:
ORF-DNA sequence BETWEEN start and stop codons. All frames would be divisible by 3.
There are 6 possible ORFs in a DNA sequence, so we have to read a sequence in 6 ways (functions in the program).

Account for incomplete strands:
has a start codon but cant find stop codon - dangling start
find a stop codon but not a start codon - dangling stop

FastAreader class - Define objects to read FastA files

OrfReader class - Turn FASTA sequence into lists of lists of Open Reading Frames (ORFs)

Methods:
init - Initialize a list of lists of ORFs that contains the frame, start, stop, and length of an ORF.

keepOrf - Collect and store attributes of the ORF for parsing later.

findOrfs - Find ORFs on a 3'-5'("top") strand. Return it as a list of ORFs.

revOrfs - Find ORFs on the 5'-3' ('bottom") strand. Return it as a list of ORFs.

revComp - creates the bottom strand of the DNA sequence


psuedocode:
slice sequence into 3 frames and iterate through them during searches

if find start codon, put into list

if start and stop codon, ORF sequence is beginning to the list to the end of stop codon

if dangling start, the whole sequence is the ORF

if dangling stop, note just the stop codon in the seq

do same thing for reverse


"""

from sequenceAnalysis3 import FastAreader, OrfFinder
########################################################################
# CommandLine
########################################################################
class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None):
        """
        Implement a parser to interpret the command line argv string using argparse.
        """
        import argparse
        self.parser = argparse.ArgumentParser(
            description='Finds the largest ORF in a DNA sequence',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store_true', help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(100, 200, 300, 500, 1000), action='store',
                                 help='minimum Gene length', default = 100)
        self.parser.add_argument('-s', '--start', action='append', nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)


########################################################################
# Main
# Here is the main program
# 
#
########################################################################
def main(inFile = None, options = None):
    """Reads in a fasta file and outputs the ORFs frame, start, stop, and length position on a output file."""
    if inFile is None:
        myCommandLine = CommandLine()
        if myCommandLine.args.longestGene:
            fastaFile = FastAreader()
            for header, sequence in fastaFile.readFasta():
                print(header)
                orfData = OrfFinder(sequence)
                orfData.findOrfs()
                orfData.revOrfs()
                filteredList = filter(lambda orf: orf[3] > myCommandLine.args.minGene, orfData.orfs)  # Filters out the ORFS depending on the minGene arg.
                for frame, start, stop, length in sorted(filteredList, key=lambda orf: orf[3], reverse = True):  # Sorts the list of ORFs by length.
                    if frame < 0:
                        print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, stop,start, length))  # flipping some of the frames so it shows up properly
                    else:    
                        print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, start,stop, length))
    else:
        myCommandLine = CommandLine(inFile)
    print(myCommandLine.args)

if __name__ == "__main__":
    main()

