function vi = GENIE3_single(expr_matrix,output_idx,input_idx,tree_method,K,ntrees)
%Computation of tree-based weights for putative edges directed towards a 
%specified target gene.
%
%vi = GENIE3_single(expr_matrix,output_idx) learns a tree model from
%expr_matrix and assigns a weight to each edge directed from a candidate
%regulator to the target gene. 
%   * expr_matrix is a matrix containing expression values. Each line
%   corresponds to an experiment and each column corresponds to a gene. 
%   * output_idx is the (column) index of the target gene in expr_matrix.
%vi is a vector of length p, where p is the number of columns in 
%expr_matrix. vi(i) is the weight of edge directed from the i-th gene of
%expr_matrix to the target gene. vi(output_idx) is set to zero.
%
%vi = GENIE3_single(expr_matrix,output_idx,input_idx) only uses as
%candidate regulators the genes whose index (as ordered in expr_matrix) is
%in input_idx. input_idx is a vector of length <= p. vi(i) such that i is
%not in input_idx is set to zero. The default vector contains the indices
%of all the genes in expr_matrix.
%
%vi = GENIE3_single(expr_matrix,output_idx,input_idx,tree_method) specifies
%which tree precedure is used. Available methods are:
%   * 'RF' - Random Forests (Default method)
%   * 'ET' - Extra Trees
%
%vi = GENIE3_single(expr_matrix,output_idx,input_idx,tree_method,K)
%specifies the number K of randomly selected attributes at each node of one
%tree. Possible values of K:
%   * 'sqrt': K = square root of the number of candidate regulators
%   (Default value)
%   * 'all': K = number of candidate regulators
%   * any numerical value
%
%vi =
%GENIE3_single(expr_matrix,output_idx,input_idx,tree_method,K,ntrees)
%specifies the number of trees grown in the ensemble. 
%Default value: 1000.

%% Check input arguments
narginchk(2,6);

if length(size(expr_matrix)) ~= 2 || size(expr_matrix,1) < 2 || size(expr_matrix,2) < 2
    error('Input argument expr_matrix must be a two-dimensional matrix, where expr_matrix[i,j] is the expression of the j-th gene in the i-th condition.')
end

if ~isscalar(output_idx)
    error('Input argument output_idx must be one integer.')
end

if ~ismember(output_idx,1:size(expr_matrix,2))
    error('Input argument output_idx must be an integer between 1 and p, where p is the number of genes in expr_matrix.')
end
   
if nargin > 2 && sum(ismember(input_idx,1:size(expr_matrix,2))) ~= length(input_idx)
    error('Input argument input_idx must be a vector containing integers between 1 and p, where p is the number of genes in expr_matrix.')
end

if nargin > 3 && sum(strcmp(tree_method,{'RF' 'ET'})) == 0
    error('Input argument tree_method must be ''RF'' (Random Forests) or ''ET'' (Extra-Trees).')
end

if nargin > 4 && ((isa(K,'char') && sum(strcmp(K,{'sqrt' 'all'})) == 0) || (isa(K,'numeric') && K <=0))
    error('Input argument K must be ''sqrt'', ''all'' or a strictly positive integer.')
end

if nargin > 5 && (~isa(ntrees,'numeric') || (isa(ntrees,'numeric') && ntrees <=0))
    error('Input argument ntrees must be a strictly positive integer.')
end

%% Data must be in single precision when used to build a tree model
expr_matrix = single(expr_matrix);
nsamples = size(expr_matrix,1); % number of experiments
ngenes = size(expr_matrix,2); % number of genes

%% Output vector
output = expr_matrix(:,output_idx);

% Normalize output data
output = output/std(output);

%% Indices of input genes
if nargin >= 3
    input_idx = unique(input_idx);
else
    % Default: all genes are candidate regulators
    input_idx = 1:ngenes;
end

% Remove target gene from candidate regulators
input_idx = setdiff(input_idx,output_idx);
ninputs = length(input_idx);


%% Tree parameters

% Default parameters: Random Forests, K=sqrt(number of input genes),
% 1000 trees in the ensemble
if nargin < 4 || (nargin >= 4 && strcmp(tree_method,'RF'))
    ok3ensparam = init_rf();
    if nargin >= 5
        if strcmp(K,'all')
            ok3ensparam=init_rf(ninputs);
        elseif isa(K,'numeric')
            ok3ensparam=init_rf(round(K));
        end
    end
elseif nargin >= 4 && strcmp(tree_method,'ET')
    ok3ensparam = init_extra_trees();
    if nargin >= 5
        if strcmp(K,'all')
            ok3ensparam=init_extra_trees(ninputs);
        elseif isa(K,'numeric')
            ok3ensparam=init_extra_trees(round(K));
        end
    end
end

% Number of trees in the ensemble
if nargin < 6
    ok3ensparam.nbterms = 1000;
else
    ok3ensparam.nbterms = ntrees;
end


%% Learning of tree model
[tree,varimp]=rtenslearn_c(expr_matrix(:,input_idx),output,int32(1:nsamples),[],ok3ensparam,0);
% Some variable importances might be slightly negative due to some rounding
% error
varimp(varimp<0) = 0;
vi = zeros(1,ngenes);
vi(input_idx) = varimp'/nsamples;
