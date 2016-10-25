#----------------------------------------------------
#datWrite.py converts standard Stacks output to json .dat files
#saves in directory 'stack_dat'
#Emily Conklin
#10/25/2016
#----------------------------------------------------

import json
import glob
import numpy as np
import os

#class for Locus object: holds ID, consensus sequence, list of SNP objects
class Locus:
    def __init__(self, ID_in, consensus_in):
        self.ID = ID_in
        self.consensus = consensus_in
        self.SNPs = []

    def getID(self):
        return self.ID

    def setSNP(self,SNPin):
        self.SNPs.append(SNPin)

    def getSNPs(self):
        return self.SNPs

#class for SNP object: holds reference nucleotide, position, list of alternate nucleotides
#& dictionary of nucleotides in each read per sample
class SNP:
    def __init__(self, ref_in, pos_in, alts_in):
        self.ref = ref_in
        self.pos = int(pos_in)
        self.alts = alts_in
        self.counts = {}

    def getAlts(self):
        return self.alts
    
    def getPos(self):
        return self.pos
    
    def getRef(self):
        return self.ref
    
    def getCounts(self):
        return self.counts

    #count tuple format:
    #{"POP*S* : {"N1":n1, "N2":n2}, POP : {} }
    def addCount(self, sample_in, seq_in):
        try:
            nucleotide = seq_in[self.pos-1]
            
            #if sample is in the dictionary
            if sample_in in self.counts.keys():
                
                #if nucleotide is already in the dictionary, add to count
                if nucleotide in self.counts[sample_in].keys():
                    self.counts[sample_in][nucleotide]+=1
                #if not, add the nucleotide at count = 1
                else:
                    self.counts[sample_in][nucleotide] = 1
        
            #if the sample's not in there yet, add at count = 1
            else:
                self.counts[sample_in] = {nucleotide : 1}

        #if file is blank, print to screen & don't include in analysis
        except:
            print(self.pos-1,len(seq_in))

####################################

#parses Stacks catalog tags.tsv file
#gets locus ID and consensus sequence to build locus objects
def makeLoci():
    with open("stacks_zip/batch_1.catalog.tags.tsv","r") as catalogTags:
        lociDict = {}
        catTagsLines = catalogTags.read().splitlines()

        #builds locus objects, adds them to lociList
        for line in catTagsLines[1:]:
            newLine = line.split()
            myLocus = Locus(newLine[2],newLine[8])
            lociDict[newLine[2]]=myLocus

    #passes {locus ID: locus obj}
    makeSNPs(lociDict)

#takes dictionary of locus objects
#parses Stacks catalog snps.tsv file
#matches locus ID to locus object, adds SNP objects
def makeSNPs(lociDict):
    with open("stacks_zip/batch_1.catalog.snps.tsv","r") as catalogSNPS:
        catSNPSLines = catalogSNPS.read().splitlines()

        for line in catSNPSLines[1:]:
            altsList = []
            newLine = line.split()
            locusID = newLine[2]
            position = newLine[3]
            reference = newLine[6]
            altsList.append(newLine[7])
            if newLine[8]!='-':
                altsList.append(newLine[8])
                if newLine[9]!='-':
                    altsList.append(newLine[9])
        
            #builds SNP object for each SNP in locus
            mySNP = SNP(reference,position,altsList)
            lociDict[locusID].setSNP(mySNP)

    #passes {locus ID: locus obj} with SNPs added
    getCounts(lociDict)

#takes dictionary of locus objects
#parses Stacks paired data files to get read counts per locus
def getCounts(lociDict):
    fasta = glob.glob("paired/*.fa")
    
    #one fasta file per locus
    for file in fasta:
        locusID = file.split("/")[1].split(".")[0]
        
        #finds matching locus object
        applicLocusObject = lociDict[locusID]
        with open(file,"r") as countFile:
            lines = countFile.read().splitlines()
            for line in lines:
                if line.startswith(">"):
                    sampleID = line.split("|")[1]
                else:
                    #add counts for each SNP in locus
                    for SNP in applicLocusObject.getSNPs():
                        SNP.addCount(sampleID, line)

    #passes {locus ID: locus obj} with counts added
    writeDatFiles(lociDict)

######################################

#dat file format:
#[{"alts":[],"counts":{"nucleotide":1},"ref":"N","pos":1}]

#takes dictionary of locus objects
#writes one .dat file per locus in directory 'stack_dat'
def writeDatFiles(lociDict):
    os.makedirs("./stack_dat/")
    for ID in lociDict.keys():
        
        SNPlist = []
        filename = "outfile_"+ID+".mod.dat"
    
        locus = lociDict[ID]   
        for SNP in locus.getSNPs():
            alts = SNP.getAlts()
            ref = SNP.getRef()
            pos = SNP.getPos()
            counts = SNP.getCounts()
            SNPstructure = {"alts" : alts, "counts": counts, "ref": ref, "pos": pos}
            SNPlist.append(SNPstructure)
        
        #filters out empty loci
        if len(SNPlist)>0:
            #converts SNPlist dictionary to json for printing
            with open(path+filename,"w") as datFile:
                datFile.write(json.dumps(SNPlist))

####################################

def main():
    makeLoci()

main()
