# HyMSMK
Integrate multiscale module kernel for disease-gene discovery in biological networks


## Requirements
Matlab 2016 or above   


## Codes 
#runCV_HyMSMK.m: cross-validation code.  <br>
This code allows parallel execution. You can change "parfor" to "for" to cancel parallel execution  <br>

#A_HyMSMK.m: the recommended version of HyMSMK in the study. <br>   
[TableScores, ~  ] = A_HyMSMK(kernel, CPF, AdjGfG,AdjGfD,AdjDfD, DisIDset, outputAB   ) <br>
% Input:  <br>
% kernel: kernel matrix for gene-gene associations<br>  
% CPF: combiniation strategy <br>
% AdjGfG: associatins between Genes (G) and Genes (G)<br>
% AdjGfD: associatins between Genes (G) and Diseases (D) <br>
% AdjDfD: associatins between Diseases (D) and Diseases (G) <br>
% DisIDset: disease index <br>  
% outputAB: output results of HN_A and HN_B <br>
% Ouput: <br>
% TableScores: a table whos variable record the scores of genes <br> 


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


 
