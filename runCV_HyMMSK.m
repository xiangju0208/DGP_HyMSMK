%     dataname = 'Dataset_HumanPPI_MencheScience2015_OMIM_HPOSim_meshDis_SymptomSim_DisGene' ;
% location for data & Results' 
datadir    = 'data';    
ResPerfDir = 'results'; if ~exist(ResPerfDir,'dir');mkdir(ResPerfDir); end
%% load  DG,GG
Com_Methodname  ='MO'   
netname         ='PPI'  
DisSetName      ='Mesh'   
fdataset        = [datadir,filesep,'demoDataSet&PPICOM_ModCM_delta=0.2.mat']
load(fdataset) 
Matrix_gene_dis00      =  AdjGfD ; 
n_disgenes_eachdisease = sum(Matrix_gene_dis00,1)';    

%% load  Kernel 
knn           = 1000      
fkernel       = [datadir,filesep,  ['k',num2str(knn),'_MO_MSK_ICCS.mat'] ]  
Kmat          = load(fkernel);   

%% 
SparseKernel  = sparse( Kmat.SparseKernel ) ; 
AdjGfG             = sparse(  AdjGfG  ); 
Matrix_gene_dis00  = sparse(  Matrix_gene_dis00  ); 
AdjDfD             = sparse(  AdjDfD  );         

% % % % % % % %  
rng(1111)          
nCVTimes  = 5; n_fold    = 5; CVtype = [num2str(n_fold),'FCV']; MinSizeDisGeSet = n_fold;    
%
dis_IDset = find(n_disgenes_eachdisease>=MinSizeDisGeSet); 
% dis_IDset = dis_IDset(1:2)   %%%for test only %%%%%%%%%%%    
n_dis =  length( dis_IDset )  

n_disease_in_Table   = length( dis_IDset   ); 
nCV_list             = zeros( n_disease_in_Table, 1 );   
matAUROC_nCVTimes    = cell(nCVTimes,1);
matAUPRC_nCVTimes    = cell(nCVTimes,1);
matRec5_nCVTimes     = cell(nCVTimes,1);
matPrec5_nCVTimes    = cell(nCVTimes,1);
methodset_nCVTimes   = cell(nCVTimes,1);
ResTime_AllCV_HO     = zeros( nCVTimes, 1 );    
ResTime_AllCV_HE     = zeros( nCVTimes, 1 );     
date_start           = datestr(now,'yyyy.mmm.dd-HH.MM.SS');
parfor  i_cv = 1:nCVTimes    
    disp(['i_cv-',num2str(i_cv) ]) 
    %
    matAUROC    = [];
    matAUPRC    = [];
    matRec5     = []; 
    matPrec5    = []; 
	tResTime_AllCV_HO = [] ; 
	tResTime_AllCV_HE = [] ; 
	methodset   = [] ;
    idx_res     = 0;     
    for ii_dis = 1:n_disease_in_Table
        tic
        Matrix_gene_dis_copy = Matrix_gene_dis00 ;  
        ID_dis               = dis_IDset(ii_dis);  
        disp(['i_cv-',num2str(i_cv),'; ii_dis-',num2str(ii_dis),'; ID_dis-',num2str(ID_dis)]) 
        ac_gene_dis00 = Matrix_gene_dis_copy(:,ID_dis ); 
        idx_pos       = find( ac_gene_dis00 );  n_pos = length( idx_pos); 
        idx_neg       = find( ~ac_gene_dis00 ); n_neg = length( idx_neg); 
        n_fold_real   = min(n_fold, n_pos) ;   
        ind_fold_pos  = crossvalind('Kfold', n_pos, n_fold_real ) ; 
        ind_fold_neg  = crossvalind('Kfold', n_neg, n_fold_real ) ; 
        for i_fold = 1:n_fold_real  
            % idx_pos_train = idx_pos(ind_fold_pos~=i_fold);
            idx_pos_test    = idx_pos(ind_fold_pos==i_fold);  n_pos_test =length(idx_pos_test); 
            %   
            idx_neg_test_WG        = idx_neg ;  n_neg_test_all = length(  idx_neg_test_WG  ) ;               
            idx_test_pos_neg_WG    = [idx_neg_test_WG; idx_pos_test ] ;  
            AdjGfD                       = Matrix_gene_dis_copy; 
            AdjGfD(idx_pos_test,ID_dis ) = 0 ;    

            %%%%%%%%%%%%% 
            TableScores1 = table; TableScores2= table;    
			% A_HyMSMK(kernel, CPF, AdjGfG,AdjGfD,AdjDfD, DisIDset, plus   ) 
			% use pre-score & Homo 
			t1=toc; %% for time estimate 			 
			[TableScores1 ] = A_HyMSMK(SparseKernel, [], AdjGfG,AdjGfD,[],     ID_dis, [1]  ) ;
            TableScores1.Properties.VariableNames = strcat(TableScores1.Properties.VariableNames,'_HO') ; 
			t2=toc; tResTime_AllCV_HO(length(tResTime_AllCV_HO)+1,1) = t2-t1;   
			
			% use pre-score & Heter   
			t1=toc;  %% for time estimate 			
			[TableScores2 ] = A_HyMSMK(SparseKernel, [], AdjGfG,AdjGfD,AdjDfD, ID_dis, [1]  ) ;  
            TableScores2.Properties.VariableNames = strcat(TableScores2.Properties.VariableNames,'_HE') ; 
			t2=toc; tResTime_AllCV_HE(length(tResTime_AllCV_HE)+1,1) = t2-t1;   

            %%%%%%%%%%%%%                              
            TableScores = [TableScores1,TableScores2]; 
            methodset = TableScores.Properties.VariableNames ;  
            %
            test_real = ac_gene_dis00(idx_test_pos_neg_WG);  
            [AUROCset, AUPRCset , Rec5set  , Prec5set ,n_method, methodnames] = getPerf(test_real,TableScores(idx_test_pos_neg_WG,:) ) ; 
            idx_res = idx_res +1 ; 
            matAUROC(idx_res,:) = AUROCset;
            matAUPRC(idx_res,:) = AUPRCset;
            matRec5(idx_res,:)  = Rec5set; 
            matPrec5(idx_res,:) = Prec5set;  
        end 
    end 
    %
	matAUROC_nCVTimes{i_cv} = mean(matAUROC,1);
	matAUPRC_nCVTimes{i_cv} = mean(matAUPRC,1);
	matRec5_nCVTimes{i_cv}  = mean(matRec5,1); 
	matPrec5_nCVTimes{i_cv} = mean(matPrec5,1); 
    methodset_nCVTimes{i_cv}= methodset;  
	ResTime_AllCV_HO(i_cv,1) = mean(tResTime_AllCV_HO,1); 
	ResTime_AllCV_HE(i_cv,1) = mean(tResTime_AllCV_HE,1);   
	toc 
    disp('  ')
% 
end 
%  
matAUROC_nCVTimes = cat(1,matAUROC_nCVTimes{:});    
matAUPRC_nCVTimes = cat(1,matAUPRC_nCVTimes{:});    
matRec5_nCVTimes  = cat(1,matRec5_nCVTimes{:}); 
matPrec5_nCVTimes = cat(1,matPrec5_nCVTimes{:}); 
% 
matRESmean        = [mean(matAUROC_nCVTimes,1);mean(matAUPRC_nCVTimes,1);mean(matRec5_nCVTimes,1);mean(matPrec5_nCVTimes,1) ]   
tbRESmean         = array2table(matRESmean, 'VariableNames', methodset_nCVTimes{1}, 'RowNames',{'AUROC','AUPRC','Rec','Prec'}) 
ResTime_AllCV_HO  = mean(ResTime_AllCV_HO,1);
ResTime_AllCV_HE  = mean(ResTime_AllCV_HE,1); 
% 
%  save 
dir_results = 'results';  
if ~exist(dir_results,'dir'); mkdir(dir_results);end 
date_cmplt  = datestr(now,'yyyy.mmm.dd-HH.MM.SS');
parastr     = sprintf('CVtype=%s_CVtime=%d_MSDGS%d', CVtype  ,  nCVTimes, MinSizeDisGeSet );   
outfile     = [dir_results,filesep,'ResPerf_HyMSMK_',Com_Methodname_abbr,'_',netname,'_',DisSetName,'_',parastr,'_',date_cmplt,'.mat'] 
save([outfile],  'tbRESmean', 'ResTime_AllCV_HO', 'ResTime_AllCV_HE', 'date_start',  'date_cmplt' ,'fkernel' , '-v7.3' )   ;   
  
