function [vecs,vals,modes,epsilon] = diffusion_map(dataMat,varargin)
% =========================================================================
% This function computes diffusion maps Coifman & Lafon 2005/2006
% On Input: dataMat  =  matrix whose ROWS correspond to observations
%           alpha    =  parameter for approximating heat kernel
%           knn      =  number of near neighbors to use to form epsilon &
%                       sparse affinity matrix
%           num_eigs =  number of eigenvectors/values to compute
%
% On output: vecs    =  eigenvectors
%            vals    =  diagonal eigenvalue matrix of T^(1/eps)
%            modes   =  coeffs for constructing data with ''vecs''
%            epsilon =  data-dependent scaling parameter for kernel:
%                               exp( -|.|^2  /  epsilon^2 )
% =========================================================================


%% input parsing
p = inputParser;

% required inputs (only the data matrix is required)
p.addRequired('dataMat', @ (x) ismatrix(x));

% parameter value iputs
p.addParameter('alpha',1/2, @isnumeric);
p.addParameter('knn',10^ceil(log10(size(dataMat,1))/2), @isnumeric);
p.addParameter('kmin', 6, @isnumeric);
p.addParameter('num_eigs',20, @isnumeric);

%% now parse the inputs
p.parse(dataMat,varargin{:});
inputs = p.Results;

disp('                                                               ')
disp('                                                               ')
disp('                                                               ')
disp('                                                               ')  
disp(inputs)

%%  define the parameters
alpha = inputs.alpha;
knn = inputs.knn;
kmin = inputs.kmin;
num_eigs = inputs.num_eigs;


%%
%%% Step 1: Set up dimensions
[nRows,nCols] = size(dataMat);

disp('--- Begin Diffusion Map Algorithm ---')



%%% build affinity

disp('      ')
disp('      ')

%%% Step 7: Find k - nearest neighbors of each y_hat in R^m

disp('      ')
disp('      ')
disp('      ')
disp('--- Begin knn search ---')
tic
[knn_idx, knn_dist] = knnsearch(dataMat , dataMat ,'K', knn);
toc
disp('--- End knn search ---')

%%% Step 8: Form a sparse (N-s) x (N-s) matrix with (N-s)k nonzero entries

%%% Step 9: Set epsilon = mean of kmin distances
epsilon = mean(knn_dist(:,kmin));
fprintf('\nepsilon = %f\n', epsilon)

%%% Step 10:  Form the sparse matrix d_hat (affinity matrix)
% kernel value at nonzero entries
dI = exp(-knn_dist.^2 / (epsilon.^2));

% prepare indices to build sparse matrix
longROW = repmat( (1:nRows)' , [knn 1] );
longCOL = knn_idx(:);

% form temporary kernel into vector
longDI = dI(:);

% create sparse array that is the matrix formed by nonlinear kernel
d_hat = sparse(longROW,longCOL,longDI, nRows, nRows);

%%% Step 11: Form symmetric Matrix--sparse matrix does not preserve
%%%                                 symmetry
size(d_hat)
J = (d_hat + d_hat') / 2;

%%% Step 12: Form the diagonal normalization matrix
P = sum(J,2);

%%% Step 13: Form the density-normalized kernel matrix
K = diag(P .^ (-alpha)) * J * diag(P .^ (-alpha));

%%% Step 14: Form the diagonal normalization matrix for new kernel
Q = sum(K,2);

%%% Step 15: Form the symmetric matrix T_hat
T_hat = diag(Q .^ (-1/2)) * K * diag(Q .^ (-1/2));

if ~issymmetric(T_hat)
    disp('T-hat not quite symmetric.  Symmetrizing.');
    T_hat = (T_hat + T_hat') ./ 2;
end
%%% Step 16: Find the L+1 largest eigenvalues/vectors
tic
disp('      ')
disp('      ')
disp('      ')
disp('--- Begin Eigen Calculations  ---')
[xi, vals] = eigs(T_hat,num_eigs);
vals = diag( vals ); % lets just keep a diagonal of eigenvalues

%%% Step 18: Compute the eigenvectors of T
vecs = diag(Q .^ (-1/2)) * xi;

disp('--- Eigenvecs Computed ---')

%%% Step 21: Compute the lth DMDC projection
modes = diag(Q) * vecs;


disp('--- DMDC Algorithm completed ---')
end

