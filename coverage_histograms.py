#----------------------------------------------------
#coverage_histograms.py parses output files from bioinformatics programs Seanome and Stacks
#obtains datasets for #SNPs per locus, #reads per SNP, and #reads per locus
#outputs integrated histograms from both datasets, either per sample or overall
#Emily Conklin
#10/25/2016
#----------------------------------------------------

import json
import seaborn as sns
import glob
import numpy as np
import matplotlib.pyplot as plt

#parses json files to get Seanome per-sample coverage
#returns list of depths per file
def seanomeReadCount():
    abundanceFilePath = "coverage_json_files"
    files = glob.glob(abundanceFilePath+"/*.json")
    allSeanome = []
    
    for file in files:
        sample = file.split(".")[0]
        
        with open(file,"r") as json:
            json = json.read()
            lines = json.split("[")
            histoData = []
            
            for term in lines[3:]:
                newTerm = term.split("]")[0].split(",")
                depth = int(newTerm[0])
                abundance = int(newTerm[1])
                for x in range(abundance):
                    histoData.append(depth)

            allSeanome.append(histoData)

    return allSeanome

#parses .tags.tsv files to get Stacks per-sample coverage
#returns list of depths per file
def stackReadCount():
    files = glob.glob("stacks_zip/POP*.tags.tsv")
    allStacks = []
    
    for file in files:
        numReads = []
        lineCounter = 0
        
        with open(file,"r") as newFile:
            count = 0
            lines = newFile.read().splitlines()
            for line in lines[1:]:
                lineCounter += 1
                parts = line.split("\t")
                type = parts[6]
                if type == "consensus" and lineCounter>2:
                    numReads.append(count)
                    count = 0
                elif type == "primary" or type == "secondary":
                    count += 1
            allStacks.append(numReads)

    return allStacks


#parses .dat files to get coverage for whole sample set
#returns count lists of SNPS per locus, reads per SNP, reads per locus
def catalogReadCount(pathIn):
    datFiles = glob.glob(pathIn+"*.dat")
    
    SNPsPerLocus_master = []
    readsPerSNP_master = []
    readsPerLocus_master = []
    
    for file in datFiles:
        #each file represents a locus
        #parses json to get data per locus
        with open(file,"r") as myFile:
            dat = myFile.read()
            parsed_dat = list(json.loads(dat))
            SNPsPerLocus = 0
            readsPerLocus = 0
            
            #data per SNP
            for SNP in parsed_dat:
                counts = SNP['counts']
                
                #filter out empty files
                if len(counts)>0:
                    SNPsPerLocus += 1
                    readsPerSNP = 0
                    
                    for sample in counts.keys():
                        nucleotides = counts[sample].keys()
                        for nuc in nucleotides:
                            reads = int(counts[sample][nuc])
                            readsPerSNP += reads
                            readsPerLocus += reads

                    readsPerSNP_master.append(readsPerSNP)

            if SNPsPerLocus>0:
                SNPsPerLocus_master.append(SNPsPerLocus)
                readsPerLocus_master.append(readsPerLocus)
    
    print("# Loci: "+str(len(SNPsPerLocus_master)))
    print("# SNPs: "+str(len(readsPerSNP_master)))
    print("Median # SNPs per locus: "+str(np.median(SNPsPerLocus_master)))
    print("Mean # SNPs per locus: "+str(np.mean(SNPsPerLocus_master)))
    print("Minimum # SNPs per locus: "+str(min(SNPsPerLocus_master)))
    print("Maximum # SNPs per locus: "+str(max(SNPsPerLocus_master)))
    print("# Reads: "+str(sum(readsPerSNP_master)))
    print("Mean # reads per SNP: "+str(np.mean(readsPerSNP_master)))
    print("Median # reads per SNP: "+str(np.median(readsPerSNP_master)))
    print("Minimum # reads per SNP: "+str(min(readsPerSNP_master)))
    print("Maximum # reads per SNP: "+str(max(readsPerSNP_master)))
    print("Mean # reads per locus: "+str(np.mean(readsPerLocus_master)))
    print("Median # reads per locus: "+str(np.median(readsPerLocus_master)))
    print("Minimum # reads per locus: "+str(min(readsPerLocus_master)))
    print("Maximum # reads per locus: "+str(max(readsPerLocus_master)))

    return SNPsPerLocus_master, readsPerSNP_master, readsPerLocus_master

#generates histogram from two datasets (superimposed)
#takes in data as lists, graph title, bin size, x-axis maximum
#saves histogram as .png file
def writeHistograms(data1,data2,displayName,binSize,maxIn):
    data1 = [float(i) for i in data1]
    data2 = [float(i) for i in data2]
    
    #generates regularly spaced bins
    bins = np.arange(0, maxIn, binSize)
    
    #writes histograms
    fig, ax = plt.subplots()
    for a in [data1, data2]:
        sns.distplot(a, bins=bins, ax=ax, kde=False)
    ax.set_xlim([0, maxIn])
    ax.set_xlabel('Read Depth')
    ax.set_ylabel('Frequency')
    ax.set_title(displayName)
    plt.savefig(displayName+'.png')
    plt.close()

#produces histograms for all samples from Seanome and Stacks
def sampleCounts():
    
    #gets count lists from both Seanome and Stacks samples
    seanome = seanomeReadCount()
    stacks = stackReadCount()

    for num in range(len(seanome)):
        if num>5:
            population = "Population 2"
        else:
            population = "Population 1"
        writeHistograms(seanome[num],stacks[num],population+" Sample "+str(num+1), 1, 20)

#produces overall histograms for Seanome and Stacks
def overallCounts():
    #seanome
    print("Seanome:")
    SNPsPerLocus_seanome, readsPerSNP_seanome, readsPerLocus_seanome = catalogReadCount("outputs/trimmed_mod_vcfs/")
    #stacks
    print("Stacks:")
    SNPsPerLocus_stacks, readsPerSNP_stacks, readsPerLocus_stacks = catalogReadCount("stack_dat/")
    
    writeHistograms(SNPsPerLocus_seanome, SNPsPerLocus_stacks, "SNPs Per Locus", 1, 20)
    writeHistograms(readsPerSNP_seanome, readsPerSNP_stacks, "Read Depth Per SNP", 5, 200)
    writeHistograms(readsPerLocus_seanome, readsPerLocus_stacks, "Read Depth Per Locus", 10, 400)

#overall histograms - comment on or off
overallCounts()
#per-sample histograms - comment on or off
sampleCounts()
