#!/usr/bin/env python3
# Name: Justin Jang (jjang12)
# Group Members: Aster Lathbury(mlathbur), Kimberly Magpantay(klmagpan), Faiz Khan(faahkhan)

"""
Read a file of FastA sequences from STDIN. Find the unique subsequences that occur in each single tRNA such that no members of this set occur among any other tRNA sets.
No member of a unique subsequence set is a substring of any other member of this set.

Example:
ACG and AAACGA are in the Unique set.

Since ACG occurs in AACGA, AACGA would be removed. ACG is Essential.

Output:
L1: tRNA name
L2: tRNA sequence with alignment characters removed
L3-onwards: each unique element

partial output example (it's way too big):
 tRNA | Glu | âˆƒUC | Bos taurus | mitochondrial
GUUCUUGU"GUUGAAUGACAACLAPGGUUUÂˆƑUCAUAPCAUUAGU?AUGGUPAG"UUCCAUGUAAGAAUACCA
GUUCUU
.UUCUUG
...CUUGU
......GU"GUU
........"GUUG
..........UUGAA
............GAAUG
..............AUGAC
................GACAAC
....................ACL
......................LA
.......................APGG
........................PGGU
...........................UUUÂ
................................ƑUC
"""

import sys
class FastAreader :
    '''Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)'''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        if self.fname == '':
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
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

class findUnique:
    ''''''
    #TA notes:
    # 3 classes needed minimum
    #powerSet, unique, essential
    def __init__(self):
        self.seqList = []  # contains list
        self.uniqueList = []  # tRNA substring list
        self.headSeqDict = {}  # sequence saving dictionary
        count = 0
        fastAFile = FastAreader()

        for head, seq, in fastAFile.readFasta():
                cleanedSequence = (seq.replace('-','').replace('.','').replace('_',''))  # clean the sequence so only sequence characters are used in the powerset
                #self.seqList.upper().strip()
    ##seqList = (CleanedSeq1, CleanedSeq2, CleanedSeq3....)
                self.headSeqDict[count] = [head, cleanedSequence]
                mypowerSet = self.powerSet(cleanedSequence)
                self.seqList.append(mypowerSet)  # After finding powerset for each tRNA, add them to this list
                count += 1

    def powerSet(self,seq):
        '''Return tRNA substrings'''
        powerSet = set()
        for i in range(len(seq)):
            size = len(seq)
            while size > i:
                powerSet.add(seq[i:size])
                size -= 1
        return powerSet  # of tRNA substrings

    def unique(self):
        '''
        Search for duplicate tRNA sequences. 
        Remove to get all unique subsets of tRNA sequences.
        '''
        for allSet in self.seqList:
            union = set()  #store union of other sets
            copySet = set()
            for powerSet in self.seqList:  # add all items from other sets to the union set
                if powerSet != allSet:  # remove common elements from the current set
                    for item in powerSet:
                        union.add(item)
            copySet = allSet - union  # Remove common elements from the current set, creating a copy

            newSet = set()
            for element in copySet:  # Add each element from the copy set to the new set
                    newSet.add(element) 
            for string1 in copySet:  # Iterate over each string in the copy set to check for possible substrings
            
                uniqueSet = set()
                for element in copySet:
                    if element != string1:
                        uniqueSet.add(element)  # Create a set of unique strings by excluding the current string
                
                for string2 in uniqueSet:  # Iterate over each unique string to compare with the current string
                    if string1 in string2:  # Check if the current string is a substring of the unique string
                        if len(string1) < len(string2):  # Check if the current string is shorter than the unique string
                            newSet.discard(string2)  # If a shorter substring is found within a longer substring, discard the longer substring

            self.uniqueList.append(newSet)  # Add the unique set to the list of unique subsets
    
    def printSeq(self):
        '''
        Prints sequences to the command line
        '''
        for index in range(len(self.headSeqDict)):  # retrieve sequence info per the current index
            mySequence = self.headSeqDict[index]
            header = mySequence[0]
            sequence = mySequence[1]

            print(header)
            print(sequence)

            size = len(sequence)  # get len of sequence
            for position in range(size):
                for substring in self.uniqueList[index]:
                    substringLength = len(substring)

                    if substring == sequence[position:position + substringLength]:  # check if the current substring matches the corresponding positionin the sequence
                        output = ('.') * position + substring  # Create a string with dots ('.') to represent positions before the substring
                        print(output)  # # Print the substring with dots representing its position
########################################################################
# Main
# Here is the main program
# 
########################################################################

def main():
    ''' '''
    #findUnique.unique()
    #findUnique.printSeq()

    tRNA = findUnique()  # instance of findUnique
    tRNA.unique()  # find subsets of tRNA sequences
    tRNA.printSeq()  # print the seuqnces and their unique subsets
    

if __name__ == "__main__":
    main()