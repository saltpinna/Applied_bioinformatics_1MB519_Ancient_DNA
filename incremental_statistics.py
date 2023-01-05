import pandas as pd
import numpy as np
import math

#---Defining variables--- Change to your need
#Telomere lengths of all chromosomes
telomeres = [3000000,3000000,0,3000000,3000000,3000000,0,3000000,3000000,3000000,3000000,3000000,3000000,3000000,300000,3000000,3000000,3000000,3000000,3000000,0,3000000,3000000,3000000,3000000,3000000,3000000,3000000]
#Name of chromosomes in file path
names = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chrX']
#File path until names part
path = 'D:/Data/mammoth.'

'''
incrementalStatistics(file_path, names, telomeres, remove_outliers)
A function to calculate the mean value and standard deviation while dividing the data and working in chunks.
Allows to calculate the new mean value and std after removing outliers.
Parts that needs to be changed for your needs is marked with *
RETURNS: A list, index 0 is the length of the file in the number of nucleotide, index 1 is the mean value, and index 2 is the std.
EXAMPLES: incrementalStatistics('../Data/mammoth.',['chr24'],[3000000]) = [24086192, 15.205945721930114, 14.45353032871932]
'''
def incrementalStatistics(filePath, names, telomeres, remove_outliers=False):
    chunk_size = 7500000

    first = False
    n = 0 #Number of nucleotides
    m = 0 #Mean value
    sq_sum = 0 #Square sum

    meanValue = 55.486431521768495 #*Used when removing outliers, this is the previous average, can be ignored if no outliers should be removed
    stdValue = 80.14920140853054 #*Used when removing outliers, this is the previous average, can be ignored if no outliers should be removed
    for i in range(len(names)):
        if telomeres[i] > chunk_size:
            telomeres[i] -= chunk_size
            first = True
        #* Loads the coverage file, change the file data type if needed and the columns that should be used
        for chunk in pd.read_csv(filePath+names[i]+'.depth.gz', sep='\t', comment='t', header=None, usecols=[4,5], skiprows=telomeres[i], chunksize=chunk_size):
            if first == False:
                chunk = chunk.mean(axis=1) #Calculates the mean value of each input row
                for index, prelMean in chunk.items(): #Iterates through each row
                    if not prelMean >= 0: #Handles possible Nan and negative values
                        prelMean = 0
                    if remove_outliers:
                        if prelMean < meanValue + 2*stdValue: #Add / remove if statement if outliers should be removed or not
                            n += 1 #Adds one nucleotide in the counter
                            m = m + ((prelMean-m)/n) #Calculates an incremental mean value
                            sq_sum += prelMean * prelMean
                    
            else:
                first = False
        print(f'Chromosome {i+1} of {len(telomeres)}')
    s = math.sqrt((sq_sum / n - m*m)) #Calculates the final std
    return [n,m,s]

print(incrementalStatistics(path,names,telomeres))