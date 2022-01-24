# HyMSMK--Integrate multiscale module kernel for disease-gene discovery in biological networks.
Biomedical data mining is very important for the research of complex diseases, and disease-gene discovery is one of the most representative topics in this field. Multiscale module structure (MMS) that widely exists in biological networks can provide useful insight for disease research. However, how to effectively mine information in MMS to enhance the ability of disease-gene discovery is challenging. Thus, we propose a type of novel hybrid methods (HyMSMK) for disease-gene discovery by integrating multiscale module kernel (MSMK) derived from multiscale module profile (MSMP). We extract MSMP with local to global structural information from comprehensive human protein interactome by multiscale modularity optimization with exponential sampling, and construct MSMK by using the MSMP as a feature matrix, combining with the relative information content of features and kernel sparsification. Then, we present several fusion strategies integrating MSMK, including a probabilistic model for rank aggregation. By a series of experiments, we study the effect of the fusion strategies and kernel sparsification on HyMSMK, and demonstrate that HyMSMK outperforms the state-of-art network-based algorithms. These results confirm that MSMK is particularly helpful for disease-gene discovery, and the kernel sparsification can improve HyMSMK in storage space and computing speed. This may provide useful insights for the study and application of MMS.      


## Requirements
Matlab 2016 or above.   


## Codes 
#runCV_HyMSMK.m: cross-validation code.  <br>
This code allows parallel execution. You can change "parfor" to "for" to cancel parallel execution  <br>
 
#A_HyMSMK.m: the recommended version of HyMSMK in the study. <br>   
[TableScores, ~  ] = A_HyMSMK(kernel, CPF, AdjGfG,AdjGfD,AdjDfD, DisIDset, outputAB   )<br>
% Input:  <br>
% kernel: kernel matrix for gene-gene associations  <br>
% CPF: combiniation strategy   <br>
% AdjGfG: associatins between Genes (G) and Genes (G)  <br>
% AdjGfD: associatins between Genes (G) and Diseases (D)   <br>
% AdjDfD: associatins between Diseases (D) and Diseases (G)   <br>
% DisIDset: disease index   <br>
% outputAB: output results of HN_A and HN_B  <br>
% Ouput: <br>
% TableScores: a table whos variable record the scores of genes.  <br> 


## Dataset
Data is located in the directory: ./data <br>
./demoDataSet_GG&DG&DD.mat, including disease-gene associations (DG), disease-disease associations (DD) and gene-gene associations (GG);  <br>
./k100_MO_MSK_ICCS.mat, including sparse kernel matrix (k=100) for gene-gene associations; <br>
./k1000_MO_MSK_ICCS.mat, including sparse kernel matrix (k=1000) for gene-gene associations. <br>


## Results 
The results will be automatically saved into the directory: results.  

## cite
If you use HyMSMK in your research, please cite: <br> 
Xiang, et al., Integrate multiscale module kernel for disease-gene discovery in biological networks.


## contact<br>
Email: xiang.ju@foxmail.com 
