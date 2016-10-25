#----------------------------------------------------
#alignment.py generates alignment .txt files from .sam files
#Emily Conklin
#10/25/2016
#----------------------------------------------------

import pysam
import glob
import numpy as np

#takes sam file as input, outputs array
def readSam(fileIn):
    samfile = pysam.AlignmentFile(fileIn, "r")

    iter = samfile.fetch()
    bigArray = []
    maxLen = 0
    
    for x in iter:
        name = x.query_name
        
        #gets reference sequence string
        seq = x.get_reference_sequence()
        #gets starting position
        pos = int(x.reference_start)
        #gets ref length, compares to initial length
        length = int(x.reference_length)
        if length > maxLen:
            maxLen = length

        #initialize empty array
        lineArray = np.empty([1,maxLen], dtype=str)
        lineArray.fill(' ')
        
        #fill array with properly aligned chars
        for num in range(pos,pos+length):
            lineArray[0,num] = seq[num-pos]
        bigArray.append(lineArray)

    writeFile(bigArray)

#writes array to file
def writeFile(dataIn):
    with open("newFile.txt","w") as aligned:
        for array in dataIn:
            aligned.write(array)
            aligned.write('\n')

def main():
    readSam(fileIn)

main()
