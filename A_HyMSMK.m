function [TableScores, ScoreMat  ] = A_HyMSMK(kernel, CPF, AdjGfG,AdjGfD,AdjDfD, DisIDset, outputAB   ) 
% Input % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % CPF: combination strategy 
% % kernel: kernel matrix for gene-gene associations 
% % AdjGfG: associatins from (f) genes (G) to Genes (G)  
% % AdjGfD: associatins from Diseases (D) to Genes (G) GfD
% % AdjDfD: associatins from Diseases (D) to Disease (G) 
% % DisIDset: disease id 
% % outputAB: outputAB=1, output the results of HN_A and HN_B.    
% % Ouput:
% % TableScores: a table whos variable record the scores of genes. 
% % 
% % % By Ju Xiang 
% % % Email: xiang.ju@foxmail.com, xiangju@csu.edu.cn      
  
	if isempty( CPF )      
		CPF = 'CPF2'; 
	end  
	if ~exist('outputAB','var') || isempty( outputAB )  
		outputAB = 0; 
	end  
    
    if length(DisIDset)>1  
		error('only allow one disease each time');	
	end 
    %  
    TableScores = table;
    sorttype    = 'descend'; 
	[N_gene, N_disease] = size( AdjGfD );  
	if ~isempty( kernel )    
		[~, scores_A ] = A_HyMSMK([], [],  AdjGfG,AdjGfD,AdjDfD, DisIDset, []  ) ;
		[~, scores_B ] = A_HyMSMK([], [],  kernel,AdjGfD,AdjDfD, DisIDset, []  ) ;  
		rmat  = getRankingOfScoreList(  [scores_A,scores_B] , sorttype, 'MeanRank' )./N_gene;   
		TableScores.('HyMSMK') = getScoreFromRankMat(rmat, CPF ) ;   
		if outputAB 
			TableScores.('RandWalk') = scores_A; 
			TableScores.('MSMK')     = scores_B;  
		end 
		
		return; 
	else       
		TableScores  = table ;   
		ac_gene_dis  = AdjGfD(:,DisIDset );  
		if isempty( AdjDfD )
			[Scores  ]= NetRandWalk(AdjGfG, 0.7, ac_gene_dis, [], [], [],  'ProbabilityNormalizationColumn'); 
			TableScores.Scores = Scores;   
		else 
			P0_D = eye( size(AdjDfD) ); P0_D = P0_D(:, DisIDset  );     
			[Scores ] = HNetRandWalk(AdjGfG,AdjGfD,AdjDfD, ac_gene_dis,P0_D, 0.7, 0.5,  0.5, 'col', []) ;
			TableScores.Scores = Scores; 	
		end  
		if nargout>=2
			ScoreMat = TableScores{:,:}; 
		end
		return; 
	end 
	% % % % % % % % % % % % % %       

end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Scores] = getScoreFromRankMat(rmat,RankMergeMethod)
	% rmat is a rank matrix, normalized by its maximal number of rows 
	% 0<= rmat(i,j) <=1, the samller, the more important 
	% Scores, the larger, the more important 
	rmin = min(rmat(:));
	rmax = max(rmat(:));
	if rmin<0 || ( rmin<1  && 1<rmax ); error('error, rmat should be 0<= rmat(i,j) <=1')
	elseif rmin>1 	
		rmat = rmat./rmax; warning( ['rmin>1, rmat is normalized by rmax'] ); 
	end 
	%
	switch RankMergeMethod
		case {'CPF1'  }    
			Scores = prod(1-rmat,2);  

		case {'CPF2' }    
			Scores = prod(1./rmat,2); 

		case {'CPF3' } 
			Scores = -sum(rmat,2);   
			
		otherwise; error(['There is no definition for the method.']);
	end

end


%% % % % % % % % % % % % % % % %  
function Yset =getExtendedInitialAssociationG2Dplus(AdjGfD,Sim_D2D, DisIDset , LogisticFuncParaC, IniProExtendType     ) 
    % Initial Association between gene_i and (query) disease_q  
    % Input: 
    % % AdjGfD    association between genes and all diseases
    % % Sim_D2D   similarity beween all diseases 
    % % DisIDset  query disease (set) (one or more diseases) 
    % % LogisticFuncParaC  
    % % ExtendType :  RandomWalkExtend2, 
    % Output
    % % Yset     extended Initial Associations between all genes and query disease(s)
    % % 
    % % By Ju Xiang 
    % % Email: xiang.ju@foxmail.com, xiangju@csu.edu.cn  
    % % 2019-3-15  
    % % % % % % % % % % % % % % % % % % % % % % % % %     
    if isempty( Sim_D2D )    
        Sim_D2D = speye(size( AdjGfD, 2));  
        warning('Simialrity is empty.!!!!!!!!!!!');
    end 
    if isempty (IniProExtendType)
        IniProExtendType = 'PrinceIniProExtend';    %  'RandomWalkExtend' 
    end       
    %
    if ~isempty (LogisticFuncParaC) 
        Sim_D2D = sparse( 1./( 1 +exp( LogisticFuncParaC.*Sim_D2D + log(9999)  )  ) );    
    end 
    % % % % %  
    switch lower( IniProExtendType )   
        case lower(  'RandomWalkIniProExtend2' ) 
            % get the clossness between Disease and all genes by by random walk
            % from the disease to adjecent disease and then to genes. 
            % % WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop )
            n_all_dis =size( Sim_D2D ,2); 
            idx = (0:n_all_dis-1)*n_all_dis+ [1:n_all_dis] ; 
            Sim_D2D( idx ) = 0 ;  % self loop = 0 
            % pro extend 
            WSim_D2D = getNormalizedMatrix(Sim_D2D, 'col', true ); 
            WAdjGfD  = getNormalizedMatrix(AdjGfD, 'col', false );   
            Yset       = WAdjGfD* WSim_D2D(:, DisIDset  ) ;    
            % two-stwp randomwalk     
            Y00set = getNormalizedMatrix(AdjGfD(:, DisIDset), 'col', false );  
            Yset = Yset + Y00set;     

        case lower( 'None' )
            Yset     = AdjGfD(:, DisIDset);   
            Yset = Yset./( sum(Yset,1) +eps );   
        otherwise
            error(['ExtendType is wrong: ', char( string( IniProExtendType ) )]);             
    end  
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %   

%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% get genescores from kernel
function [scores ] = getGeneScoreByKernel( kernel, ac_gene_dis, KernelScoreMethod )
% % % % % % % % % % % 
    switch KernelScoreMethod
        case 'RWR'
            scores = A_RWRplus( kernel, 0.7, ac_gene_dis, [], [],  [],  []) ;

        case 'DN' % use sum of similarity scores to seed genes. 
            scores = kernel*ac_gene_dis ; 

        case 'RLS' % use  it can be improved     
            scores = kernel*( kernel+1*eye( size(kernel) ) )^-1*ac_gene_dis ; 

        otherwise; error('There is no defintion');
    end  
 
end   


%% % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [Ranks, ord ]= getRankingOfScoreList(   ScoreList, sorttype, OP_IdenticalValues )
    % by Xiang    
    % 2019-2 
    if ~exist('OP_IdenticalValues','var') || isempty(OP_IdenticalValues)
    % for elements with the same values 
       OP_IdenticalValues = 'MeanRank';  
    end
    %
    [~, ord ] = sort(ScoreList , 1,  sorttype );  
    
    IsSparse_ScoreList = issparse( ScoreList );
    if IsSparse_ScoreList   
        ScoreList = full( ScoreList ); 
    end
    %
    IDlist = [1: size( ScoreList, 1 ) ]'; 
    Ranks = zeros( size( ScoreList ) ); 
    rank_t = zeros( size(IDlist) ); 
    for d2=1:size( ScoreList ,2 )
        for d3 = 1:size( ScoreList ,3)
            for d4 = 1:size( ScoreList ,4)
                score_t    = ScoreList(:,d2,d3,d4); 
                rank_t(ord(:,d2,d3,d4)) = IDlist; 
                % Ranks( ord(:,d2,d3,d4),   d2,d3,d4) = IDlist ; 
                % 
                if ~strcmpi(OP_IdenticalValues, 'None')    
                    [uniqueScores, ~,ic] = unique( score_t ) ;
                    if length( ic  ) ~= length( uniqueScores  )
                        for ii_uniqueScorese = 1: length( uniqueScores  ) 
                            idx = ( ic== ii_uniqueScorese ) ; 
                            n_thisscore = nnz( idx )   ;
                            if n_thisscore>1
                                if strcmpi(OP_IdenticalValues, 'MeanRank') 
                                    % Ranks( idx,   d2,d3,d4 ) = mean(Ranks(idx,   d2,d3,d4), 1);  %对于相同分值的元素赋予相同的所有ranks的均值
                                    rank_t( idx ) = mean(rank_t(idx), 1);  %对于相同分值的元素赋予相同的所有ranks的均值
                                    % % sum( labels_ord( idx ) )/nnz( idx )
                                elseif strcmpi(OP_IdenticalValues, 'RandPermutation') 
                                    ind = find(idx); 
                                    ind_randperm = ind(  randperm( n_thisscore  )  ); 
                                    rank_t( ind ) = rank_t( ind_randperm ); 
                                else
                                    error('There is no definition of OP_IdenticalValues');
                                end
                            end  
                            % %sum(labels_ord) 
                        end
                    end
                end
                Ranks( : ,   d2,d3,d4) = rank_t;  
                % % 
            end
        end
    end 
    
    if IsSparse_ScoreList   
        Ranks = sparse( Ranks );  
    end

end

%% 
function WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop ) 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% Adj  adjecent matrix
% % NormalizationType: 
% % 'probability normalization' for 'Random Walk' RWR, RWRH    
% % 'laplacian normalization' for prince and more....
% SetIsolatedNodeSelfLoop    set isolated node
% >= Matlab 2016
% % % % % % % % % % % % % % % % % % % % % % % % % 
    if ischar(NormalizationType)
    %         NormalizationType =  (NormalizationType);
        switch  lower( NormalizationType )
            case lower( { 'column','col',  ...
                    'ProbabilityNormalizationColumn','ProbabilityNormalizationCol',...
                    'ProbabilityColumnNormalization','ProbabilityColNormalization',...
                    'NormalizationColumn','NormalizationCol' , ...
                    'ColumnNormalization','ColNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =1;
            case lower({ 'row' ,'ProbabilityNormalizationRow' ,'NormalizationRow' ,'ProbabilityRowNormalization' ,'RowNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =2;
            case lower('LaplacianNormalization')
                NormalizationName = NormalizationType; 
            case lower('LaplacianNormalizationMeanDegree')
                NormalizationName = NormalizationType;  
            case lower({'none', 'None', 'NONE'})
                % NormalizationName = 'None'; 
                WAdj = Adj; 
                return; 
            otherwise
                error(['There is no type of normalization: ',char( string(NormalizationType) )] );
        end
        
    elseif isnumeric(  NormalizationType   ) 
        NormalizationName =  ( 'ProbabilityNormalization' ) ;  %  'Random Walk'  
        dim = NormalizationType; 
        
    elseif isempty( NormalizationType )
        WAdj = Adj; 
        return;  
        
    else; error('There is no defintion of NormalizationType')
    end  
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     %
    switch lower( NormalizationName )
        case lower( 'ProbabilityNormalization' )
            degrees = sum(Adj,dim);
            if any( degrees~=1)
                WAdj = Adj./ ( degrees+eps  );            
            else
                WAdj = Adj; 
            end
            % 
            if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2) 
                ii = find( ~degrees ); 
                idx = sub2ind( size(Adj), ii,ii ); 
                WAdj(idx) = 1;  % set to be 1 for isolated nodes, 
            end
            
        case lower( 'LaplacianNormalization')
            deg_rowvec = ( sum(Adj,1) ).^0.5;  
            deg_colvec = ( sum(Adj,2) ).^0.5;   
            WAdj = (Adj./(deg_colvec+eps))./(deg_rowvec+eps) ;    
            % 
            if SetIsolatedNodeSelfLoop && size(Adj,1)==size(Adj,2)
                ii = find( ~sum(Adj,2) ) ;  
                WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
            end 
            
        case lower( {'None','none'} )
            WAdj = Adj;   
            
        otherwise
            error(['NormalizationName is wrong: ',char(string(NormalizationName) )   ]);
    end
% % % % 
end


%% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function  [Pt, WAdj ]= NetRandWalk(Adj, r_restart, P0, N_max_iter, Eps_min_change, IsNormalized,  NormalizationType)    
% need getNormalizedMatrix  
% Input 
% Adj
% r_restart
% P0     it should be normalized 
% N_max_iter
% Eps_min_change
% IsNormalized   Whether Adj has been normalized: True or False  
% NormalizationType:  
% (1) random walk with restart {'ProbabilityNormalizationColumn','ProbabilityNormalizationCol','col','column'}        
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Ouput
% Pt  
% WAdj   normalized Adj 
% % % % % % % % % % % % %   
    if ~exist('N_max_iter','var') || isempty(N_max_iter) || (isnumeric( N_max_iter) && N_max_iter<=1 ) 
        N_max_iter =100; 
    elseif ~isnumeric( N_max_iter)  
        error('N_max_iter should be isnumeric!!!!!') ;
    end
    %
    if ~exist('Eps_min_change','var') || isempty(Eps_min_change) 
        Eps_min_change =10^-6; 
    elseif isnumeric( Eps_min_change) && Eps_min_change>=1 
        warning('The Eps_min_change is nomenaning. Reset Eps_min_change to be 10^-6.'); 
        Eps_min_change =10^-6;  
    elseif ~isnumeric( Eps_min_change)  
        error('Eps_min_change should be isnumeric!!!!!') ;
    end
    
    if ~exist('IsNormalized','var') || isempty(IsNormalized) 
        IsNormalized = false;  % Adj has been normalized for fast run.   
    end
    
    if ~exist('NormalizationType','var') || isempty(NormalizationType) 
        NormalizationType = 'ProbabilityNormalizationColumn'; %%for 'Random Walk' RWR, RWRH  and more   
    end        
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    P0  = full(P0); 
    %     
    if IsNormalized 
        WAdj = Adj; 
    else 
        switch NormalizationType
            case {'ProbabilityNormalizationColumn','ProbabilityNormalizationCol','col','column'} 
                WAdj = getNormalizedMatrix(Adj, 'col', true );
                %%%P0 = P0./(sum(P0,1)+eps);    % total probability is 1.    
                
            otherwise
                error(['NormalizationType is wrong: ',char( string(NormalizationType) )]); 
        end        
    end
   
    % % Solver_IterationPropagation
    % % It can be used directly when IsNormalized is TRUE.  
    Pt = P0;
    for T = 1: N_max_iter
        Pt1 = (1-r_restart)*WAdj*Pt + r_restart*P0;
        if all( sum( abs( Pt1-Pt )) < Eps_min_change )
            break;
        end
        Pt = Pt1;
    end
    Pt = full(Pt);    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
end 

  
%%  
function [P_G, P_D, Pt ] = HNetRandWalk(AdjGfG,AdjGfD,AdjDfD, P0_G,P0_D, restart, pro_jump, eta, NormalizationType, isdebug) 
% % % % % % % % % % % % % % % % % % % %   
% Input % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% AdjGfG: associatins from (f) genes (G) to Genes (G)  
% AdjGfD: associatins from Diseases (D) to Genes (G) GfD
% AdjDfD  associatins from Diseases (D) to Disease (G) 
% P0_G: column vector (set) initial probabilities in Gene network
% P0_D: column vector (set) initial probabilities in Disease network
% P0_G and P0_D must have the same # of columns. 
% gamma/restart: restarting Probability  
% pro_jump: jumping Probability between different networks
% eta: ratio of Probability in second network
% NormalizationType = 'ProbabilityNormalizationColumn';     
% Ouput % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% P_G: stable probabilities in Gene network 
% P_D: stable probabilities in Disease network  
% Pt: stable probabilities in Gene+disease network 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % NetRandWalk
% random walk with restart algorithm in heterogeneous network. 
% Reference: Yongjin Li, Jagdish C. Patra, Genome-wide inferring gene–phenotype relationship 
% by walking on the heterogeneous network, Bioinformatics, 2010, 26: 1219-1224.  
    
    if ~exist('pro_jump','var') || isempty (pro_jump)
        pro_jump = 0.5; 
    end   
    if ~exist('eta','var') || isempty (eta)
        eta = 0.5; 
    end   
    if ~exist('NormalizationType','var') || isempty (NormalizationType)
        NormalizationType = 'ProbabilityNormalizationColumn' ; 
    elseif ~ismember(NormalizationType,{  'column','col',  'ProbabilityNormalizationColumn','ProbabilityNormalizationCol'     } )
        error(['NormalizationType is wrong: ',char(string(NormalizationType)) ]);
    end   
   if ~exist('isdebug','var') || isempty (isdebug)
        isdebug = false;  
    end   
    %  
    [N_gene, N_disease] = size( AdjGfD );
    if isempty( AdjDfD )
       AdjDfD = speye(N_disease);  
       warning('AdjDfD is empty.');
    end
    % 
    if size(P0_G,1)~=N_gene; error( 'P0_G must be column vector(set), with length of the number of genes.'  );end
    if size(P0_D,1)~=N_disease; error( 'P0_D must be column vector(set), with length of the number of diseases.'  );end
    %
    [ M , IsNormalized ] = getNormalizedMatrix_Heter(AdjGfG,AdjGfD,AdjDfD, pro_jump,  NormalizationType, []) ; 
  
    if any( strcmpi(NormalizationType,{'col','ProbabilityNormalizationColumn','NormalizationColumn', 'Column'}) )
        P0_G = P0_G./(sum(P0_G,1)+eps);  
        P0_D = P0_D./(sum(P0_D,1)+eps); 
        P0 = [ (1-eta)*P0_G; eta*P0_D]; 
		P0 = P0./( sum(P0,1)+eps )  ;   
    else
        P0 = [ (1-eta)*P0_G; eta*P0_D];
    end     
    Pt = NetRandWalk(M, restart, P0 , [],[], IsNormalized);   
    P_G = Pt(1:N_gene    ,:);
    P_D = Pt(N_gene+1:end,:); 
    
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [ combMatrix , IsNormalized ] = getNormalizedMatrix_Heter(AdjGfG,AdjGfD,AdjDfD, pro_jump,  NormalizationType, isdebug) 
% % % % % % % % % % % % % % % % % % % %  
% Input % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% AdjGfG: associatins from (f) genes (G) to Genes (G)  
% AdjGfD: associatins from Diseases (D) to Genes (G) GfD
% AdjDfD  associatins from Diseases (D) to Disease (G)  
% pro_jump: jumping Probability from first layer to second layer or weighting the effect of second layer on the first layer.    
% NormalizationType = 'ProbabilityNormalizationColumn'; %%   
% NormalizationType = 'None'; %%  without normalization ....    
% Ouput % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% combMatrix is matrix after normalization.   
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% By Ju Xiang, 
% Email: xiang.ju@foxmail.com, xiangju@csu.edu.cn  
% 2019-8-2 
    % global   Global_Var_RWRH 
    if ~exist('pro_jump','var') || isempty (pro_jump)
        pro_jump = 0.5; 
    elseif pro_jump>1 || pro_jump <0
        error('pro_jump is wrong. it should be between 0 and 1');
    end       
   if ~exist('isdebug','var') || isempty (isdebug)
        isdebug = false;  
    end   
    %  
    [N_gene, N_disease] = size( AdjGfD );
    if isempty( AdjDfD )
       AdjDfD = speye(N_disease);  
       warning('AdjDfD is empty.');
    end 
    %
    if ~exist('NormalizationType','var') || isempty(NormalizationType)
        NormalizationType = 'None'; 
    end
    %  
    IsNormalized = true; 
    switch lower( NormalizationType )
        case lower( {'None'} )
            combMatrix = [ AdjGfG, AdjGfD; AdjGfD', AdjDfD    ] ;
            IsNormalized = false;
            
        case lower( {'Weight'} )  
            combMatrix = [ (1-pro_jump).*AdjGfG, pro_jump.*AdjGfD; pro_jump.*AdjGfD', (1-pro_jump).*AdjDfD    ] ;
            IsNormalized = false;
            
        case lower( {'col','ProbabilityNormalizationColumn','NormalizationColumn', 'Column'} )   %概率解释  ，确保列和为1 
            idxDis_WithDiseaseGene =  sum( AdjGfD, 1)~=0;   % mark diseases with disease-genes
            idxGene_WithDisease    = (sum( AdjGfD, 2)~=0)';   % mark genes that are associated with diseases
            % WAdj = getNormalizedMatrix(Adj, NormalizationType, SetIsolatedNodeSelfLoop )
            M_GfG = getNormalizedMatrix(AdjGfG   , NormalizationType, true  ); 
            M_DfD = getNormalizedMatrix(AdjDfD   , NormalizationType, true  ); 
            M_GfD = getNormalizedMatrix(AdjGfD   , NormalizationType, false );  % probabilities from disease space to gene space 
            M_DfG = getNormalizedMatrix(AdjGfD'  , NormalizationType, false );  % probabilities from gene space to disease space
            %
            M_GfG(:,idxGene_WithDisease)       = (1-pro_jump).*M_GfG(:,idxGene_WithDisease); 
            M_DfD(:,idxDis_WithDiseaseGene )   = (1-pro_jump).*M_DfD(:,idxDis_WithDiseaseGene ) ; 
            M_GfD                           = pro_jump.*M_GfD; % Disease-columns without disease-genes is all zeros. So no use idxDis_WithDiseaseGene
            M_DfG                           = pro_jump.*M_DfG; % Gene-columns without diseases is all zeros. So no use idxGene_WithDisease
            %
            combMatrix = [ M_GfG, M_GfD; M_DfG, M_DfD    ] ; 
                        
        otherwise
            error('No definition.');
    end
    
    %  
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
