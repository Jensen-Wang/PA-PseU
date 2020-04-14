# -*- coding: utf-8 -*-

import numpy as np

# Read physicochemical properties
def loadTXTfile():
    tmp = np.loadtxt("physichemical.csv", dtype=np.float, delimiter=",")
    return tmp

physichemical = loadTXTfile()
print(physichemical.shape)

# Perform a standard conversion
physichemical_stand = np.zeros((10,16), dtype=float)
for i in range(physichemical.shape[0]):
    for j in range(physichemical.shape[1]):
        physichemical_stand[i][j] = (physichemical[i][j] - np.mean(physichemical[i]))/np.std(physichemical[i]) 
        
# List of the 16 different dinucleotides
list = ["AA","AC","AG","AU","CA","CC","CG","CU","GA","GC","GG","GU","UA","UC","UG","UU"]
# Calculate the functiion of Lambda for every RNA sequence
def getH(seq, j, physichemical_stand,i):
    dinu1 = seq[i:i+2]
    dinu1_p = list.index(dinu1)
    dinu2 = seq[i+j:i+j+2]
    dinu2_p = list.index(dinu2)
    sum_u = 0
    for u in range(physichemical_stand.shape[0]):
        dinu1_physi_value = physichemical_stand[u][dinu1_p]
        dinu2_physi_value = physichemical_stand[u][dinu2_p]
        sum_u = sum_u + pow(dinu1_physi_value-dinu2_physi_value,2)
    return sum_u / (physichemical_stand.shape[0])

# Get the sequence order correlated factors
def getThetaLambda(seq, Lambda, physichemical_stand):
    sum_theta = 0
    for i in range(len(seq)-1-Lambda):
        sum_theta = sum_theta + getH(seq, Lambda, physichemical_stand, i)
    return sum_theta/(len(seq)-1-Lambda)
# Store all the RNA sequences
all_RNA_sequences = []

##The code need to be ran by batch for getting different file for every data set
## Read the RNA sequences and add these RNA sequences into all_RNA_sequences 
#with open("SuppS1.txt", "r") as file_object:
#    S1 = file_object.readlines()
#    all_RNA_sequences.extend(S1)
#with open("SuppS2.txt", "r") as file_object:
#    S2 = file_object.readlines()
#    all_RNA_sequences.extend(S2)
#with open("SuppS3.txt", "r") as file_object:
#    S3 = file_object.readlines()
#    all_RNA_sequences.extend(S3)
#with open("SuppS4.txt", "r") as file_object:
#    S4 = file_object.readlines()
#    all_RNA_sequences.extend(S4)
with open("SuppS5.txt", "r") as file_object:
    S5 = file_object.readlines()
    all_RNA_sequences.extend(S5)
# Calculate every theta
def getThetaList(all_RNA_sequences,lambda_scale,physichemical_stand):
    num = 0
    theta_list = []
    for seq in all_RNA_sequences:
        if seq.startswith(">"):
            num = num + 1
        else:
            seq = seq.rstrip()
            theta_list_line = []
            for Lambda in range(1,lambda_scale+1):
                theta_list_line.append(getThetaLambda(seq, Lambda, physichemical_stand))
            theta_list.append(theta_list_line)
    return np.asarray(theta_list)

# Get sequence order correlated factors for lambda
theta_list = getThetaList(all_RNA_sequences,29,physichemical_stand) #Here, 29 for S2, S5 and 19 for S1, S3, S4. It is changed according to length of RNA sequence
print(theta_list.shape)
#theta_list = getThetaList(all_RNA_sequences,19,physichemical_stand)
#print(theta_list.shape)

# Store the sequence order correlated factors into data_theta.npy
#theta_list = np.asarray(theta_list)
#np.save("HSapiensSeqOrder.npy", theta_list)

theta_list = np.asarray(theta_list)
np.save("IndeptSCerevisiaeSeqOrder.npy", theta_list)

#theta_list = np.asarray(theta_list)
#np.save("IndptHSapiensSeqOrder.npy", theta_list)