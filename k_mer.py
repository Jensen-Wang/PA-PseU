# -*- coding: utf-8 -*-


import itertools

bases=["A", "U", "G", "C"]
def create_features(k):
    list = ["".join(p) for p in itertools.product(bases, repeat=k)]
    return list

def create_dict(list):
    kmers = {}
    for i in list:
        kmers[i]=0
    return kmers
import numpy as np



# Get the k-mer nucleotide composition
#the read code need to be implemented in batches by S1,S2,S3 to generate different 8-mer files for different species
def getArray(k):
    # 944 is the total number of RNA sequences in S3
    # the number should be changeed to 628 for total number of RNA sequences in S2
    # the number should be changeed to 990 for total number of RNA sequences in S1
    
    k_array = np.zeros((944, pow(4,k)), dtype=int)
    # Get 4^k different k-mer nucleotide composition
    kmers = create_dict(create_features(k))
    # Read the RNA sequences
#    with open("C:/Users/Jensen Wang/Desktop/Supp/SuppS1.txt", "r") as file_object:
#        SuppS1 = file_object.readlines()
#    with open("C:/Users/Jensen Wang/Desktop/Supp/SuppS2.txt", "r") as file_object:
#        SuppS2 = file_object.readlines()
    with open("C:/Users/Jensen Wang/Desktop/Supp/SuppS3.txt", "r") as file_object:
        SuppS3 = file_object.readlines()
    
      
    # Get the occurrence number of each k-mer for a RNA sequence
    integer = -1
    for sequences in [SuppS3]: #Here, variable should be chaanged for different RNA of species
        num = 0
        for seq in sequences:
            for i in kmers:
                kmers[i] = 0
            num += 1
            if num%2 == 0:
                integer += 1
                seq = seq.rstrip()
                for i in range(len(seq) - k + 1):
                    kmer = seq[i:i+k]
                    if kmer in kmers.keys():
                        kmers[kmer] += 1
                j = -1
                for kmer, count in kmers.items():
                    j += 1
                    k_array[integer][j] = count
   
    
    
     
    # Calculate the occurrence frequency of each k-mer for a RNA sequence
    k_array = k_array.astype("float")
    for i in range(k_array.shape[0]):
        k_array[i] = k_array[i] / sum(k_array[i])
    return k_array





def getIndependentArray(k):
    # 200 is the number of RNA sequences for Independent data sets S4 or S5
    k_array = np.zeros((200, pow(4,k)), dtype=int)
    # Get 4^k different k-mer nucleotide composition for Independent data sets
    # data sets of S4,S4 name should be changed by batches
    kmers = create_dict(create_features(k))
    # Read the RNA sequences
#    with open("C:/Users/Jensen Wang/Desktop/Supp/SuppS4.txt", "r") as file_object:
#        SuppS4 = file_object.readlines()
    with open("C:/Users/Jensen Wang/Desktop/Supp/SuppS5.txt", "r") as file_object:
        SuppS5 = file_object.readlines()
        
    integer = -1
    for sequences in [SuppS5]: #SuppS5 should be changed to SuppS4 when calculating S4
        num = 0
        for seq in sequences:
            for i in kmers:
                kmers[i] = 0
            num += 1
            if num%2 == 0:
                integer += 1
                seq = seq.rstrip()
                for i in range(len(seq) - k + 1):
                    kmer = seq[i:i+k]
                    if kmer in kmers.keys():
                        kmers[kmer] += 1
                j = -1
                for kmer, count in kmers.items():
                    j += 1
                    k_array[integer][j] = count
     
    # Calculate the occurrence frequency of each k-mer for a RNA sequence
    k_array = k_array.astype("float")
    for i in range(k_array.shape[0]):
        k_array[i] = k_array[i] / sum(k_array[i])
    return k_array





# Set the appropriate label for each class of RNA sequences
def getTarget():
    # RNA sequences 472-944 in S3 contains no pseudouridine sites, and we set the label of these sequences as -1
    # RNA sequences 314-628 in S2 contains no pseudouridine sites, and we set the label of these sequences as -1
    # RNA sequences 495-990 in S1 contains no pseudouridine sites, and we set the label of these sequences as -1
    
    target = np.ones((944,), dtype=int)
#    for i in range(495,990): 
#        target[i] = -1

#    for i in range(314, 628):
#        target[i] = -1
#
    for i in range(472, 944):
        target[i] = -1
    return target

def getIndependentTarget():
    # RNA sequences 100-200 in S4 contains no pseudouridine sites, and we set the label of these sequences as -1
    # RNA sequences 100-200 in S5 contains no pseudouridine sites, and we set the label of these sequences as -1
    
    target = np.ones((200,), dtype=int)
    for i in range(100,200): 
        target[i] = -1



    return target





# The following file generator code also need to be implemented by batches matching with above code.


#human8mer = getArray(8)
#print(human8mer.shape)
#human8mer = np.asarray(human8mer)
#np.save("HSapiens8mer.npy", human8mer)


#SCerevisiae8mer = getArray(8)
#print(SCerevisiae8mer.shape)
#SCerevisiae8mer = np.asarray(SCerevisiae8mer)
#np.save("SCerevisiae8mer.npy", SCerevisiae8mer)

#MMusculus8mer = getArray(8)
#print(MMusculus8mer.shape)
#MMusculus8mer = np.asarray(MMusculus8mer)
#np.save("MMusculus8mer.npy", MMusculus8mer)

#IndependentHSapiens8mer = getIndependentArray(8)
#print(IndependentHSapiens8mer.shape)
#IndependentHSapiens8mer = np.asarray(IndependentHSapiens8mer)
#np.save("IndependentHSapiens8mer.npy", IndependentHSapiens8mer)

#IndependentSCerevisiae8mer = getIndependentArray(8)
#print(IndependentSCerevisiae8mer.shape)
#IndependentSCerevisiae8mer = np.asarray(IndependentSCerevisiae8mer)
#np.save("IndependentSCerevisiae8mer.npy", IndependentSCerevisiae8mer)
    
##Store the labels as target and IndependentTarget
target = np.asarray(getTarget())
np.save("TargetMMusculus.npy", target)

#target = np.asarray(getTarget())
#np.save("TargetHSapiens.npy", target)

#target = np.asarray(getTarget())
#np.save("TargetSCerevisiae.npy", target)


IndependentTarget = np.asarray(getIndependentTarget())
np.save("IndeptTargetSCerevisiae.npy",IndependentTarget)


#IndependentTarget = np.asarray(getIndependentTarget())
#np.save("IndeptTargetHSapiens.npy",IndependentTarget)
