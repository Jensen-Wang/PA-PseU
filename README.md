# PA-PseU
## This repo is for the data and code for my paper "PA-PseU: An Incremental Passive-Aggressive Based Method For Identifying RNA Pseudouridine Sites".  
----------------------------------------------------
Pseudouridine (__Ψ__) is the most abundant RNA modification existing in ubiquitous organisms and participates a bunch of biological processes. Identification of pseudouridine sites have significant meanings for the study of biological process and drug development related to RNA. Wet experiments for pseudouridine sites are expensive, time-consuming and easily influenced by the environment although they can detect __Ψ__ sites in whole transcriptome. Although some of the computational methods have been presented, their performance are still unsatisfactory and their computation cost are usually large. In this article, we propose an incremental identification method called PA-PseU which based on Passive-Aggressive algorithm as __Ψ__ sites classifier. The combination of chi-square test and logistic regression were used to select the optimal feature subsets. Its effectiveness has been demonstrated via 10-fold cross-validations, jack-knife tests and independent tests. Results demonstrated that our model outperformed the state-of-art methods by a big margin and improved the previous models substantially. Moreover, PA-PseU is a high-efficient method which take low computation cost with rapid computation speed. We provide a paradigm of constructing a well performance model which can be applied in large amount of data when other methods do not work. It can be anticipated that PA-PseU will be an useful tool for identifying the __Ψ__ sites in genome analysis.

**`SuppS1.txt`** is the benchmark dataset for H. sapiens.  

**`SuppS2.txt`** is the benchmark dataset for S. cerevisiae.  

**`SuppS3.txt`** is the benchmark dataset for M. musculus.  

**`SuppS4.txt`** is the independent dataset for H. sapiens.  

**`SuppS5.txt`** is the benchmark dataset for S. cerevisiae.  

**`k_mer.py`** is for the k-mer feature extraction. You can get the k-mer feature of any sequence for the training and testing datasets by running the script.  

**`SequenceOder.py`** is for the sequence feature extraction. You can get the sequence feature of any sequence for the training and testing datasets by running the script.

**`model.ipynb`** is the main code of the mdoel which can be run in jypyter notebook. It contains the results and pictures of our model. 



