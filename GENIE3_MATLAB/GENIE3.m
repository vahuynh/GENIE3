function VIM = GENIE3(expr_matrix,input_idx,tree_method,K,ntrees)
%Computation of tree-based weights for all putative edges.
%
%VIM = GENIE3(expr_matrix) learns p tree models from expr_matrix, where p 
%is the number of columns (genes) in expr_matrix, and assigns a weight to 
%each edge directed from any gene in expr_matrix to any other gene. 
%expr_matrix is a matrix containing expression values. Each line 
%corresponds to an experiment and each column corresponds to a gene. VIM is
%a matrix of size p x p. VIM(i,j) is the weight of edge directed from the
%i-th gene of expr_matrix to the j-th gene. VIM(i,i) is set to zero for all
%i.
%
%VIM = GENIE3(expr_matrix,input_idx) only uses as candidate regulators the
%genes whose index (as ordered in expr_matrix) is in input_idx. input_idx
%is a vector of length <= p. VIM(i,:) such that i is not in input_idx is
%set to zero. The default vector contains the indices of all the genes in
%expr_matrix.
%
%VIM = GENIE3(expr_matrix,input_idx,tree_method) specifies
%which tree precedure is used. Available methods are:
%   * 'RF' - Random Forests (Default method)
%   * 'ET' - Extra Trees
%
%VIM = GENIE3(expr_matrix,input_idx,tree_method,K)
%specifies the number K of randomly selected attributes at each node of one
%tree. Possible values of K:
%   * 'sqrt': K = square root of the number of candidate regulators
%   (Default value)
%   * 'all': K = number of candidate regulators
%   * any numerical value
%
%VIM =
%GENIE3(expr_matrix,input_idx,tree_method,K,ntrees)
%specifies the number of trees grown in each ensemble. 
%Default value: 1000.


%%
tic;

%% Check input arguments
narginchk(1,5);

if length(size(expr_matrix)) ~= 2 || size(expr_matrix,1) < 2 || size(expr_matrix,2) < 2
    error('Input argument expr_matrix must be a two-dimensional matrix, where expr_matrix(i,j) is the expression of the j-th gene in the i-th condition.')
end

ngenes = size(expr_matrix,2);

if nargin > 1 && sum(ismember(input_idx,1:ngenes)) ~= length(input_idx)
    error('Input argument input_idx must be a vector containing integers between 1 and p, where p is the number of genes in expr_matrix.')
end

if nargin > 2 && sum(strcmp(tree_method,{'RF' 'ET'})) == 0
    error('Input argument tree_method must be ''RF'' (Random Forests) or ''ET'' (Extra-Trees).')
end

if nargin > 3 && ((isa(K,'char') && sum(strcmp(K,{'sqrt' 'all'})) == 0) || (isa(K,'numeric') && K <=0))
    error('Input argument K must be ''sqrt'', ''all'' or a strictly positive integer.')
end

if nargin > 4 && (~isa(ntrees,'numeric') || (isa(ntrees,'numeric') && ntrees <=0))
    error('Input argument ntrees must be a strictly positive integer.')
end


%% Print parameters
if nargin < 3
    print_method = 'RF';
else
    print_method = tree_method;
end
    
if nargin < 4
    print_K = 'sqrt';
else
    print_K = K;
    if isa(K,'numeric')
        print_K = num2str(K);
    end
end

if nargin < 5
    print_ntrees = 1000;
else
    print_ntrees = ntrees;
end

fprintf('Tree method: %s\n',print_method)
fprintf('K: %s\n',print_K)
fprintf('Number of trees: %d\n\n',print_ntrees)


%% Learn an ensemble of trees for each gene and compute importance


VIM = zeros(ngenes,ngenes);

for i=1:ngenes
    
    fprintf('Gene %d/%d...\n',i,ngenes)
    
    if nargin == 1
        VIM(i,:) = GENIE3_single(expr_matrix,i);
    elseif nargin == 2
        VIM(i,:) = GENIE3_single(expr_matrix,i,input_idx);
    elseif nargin == 3
        VIM(i,:) = GENIE3_single(expr_matrix,i,input_idx,tree_method);
    elseif nargin == 4
        VIM(i,:) = GENIE3_single(expr_matrix,i,input_idx,tree_method,K);
    else
        VIM(i,:) = GENIE3_single(expr_matrix,i,input_idx,tree_method,K,ntrees);
    end
end

VIM = VIM';

%% 
toc;