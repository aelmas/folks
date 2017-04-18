function [F,models] = mappingFunction(T,M,D,verbose)
% mappingFunction Computes the enrichment of the sequence features relative to
%   their background distribution.
%   [F,models] = mappingFunction(T,M,D,verbose) Computes the enrichment of the 
%   sequence model(s) in the given sequence set T relative to the 
%   enrichments in the background (D).
%
%   T, a set of character strings representing DNA sequences, or is the 
%     path to the fasta file containing DNA sequences.
%
%   M, a matrix where each row represents a gapped n-mer model.
%   Example
%     [1 1 1] represents all 3-mers ['AAA','AAC',...,'ACA','ACC',...,'TTT']
%     [1 0 1] represents the gapped 3-mers where the 2nd nucleotide
%       is ignored ['AgA','AgC',...,'CgA','CgC',...,'TgT']
%     Note: [1 0 1] have the same number of elements with 2-mers.
%
%   D, struct array. [model] Each row contains the corresponding sequence features, 
%     e.g., (gapped) n-mer distributions. [distribution] Contains the
%     corresponding distribution of the model.
%
%   verbose, logical value for displaying the progress. default is false.
%
%   F, enrichment vectors representing the relative enrichment of the given 
%     sequence features in rows.
%   
%   models, all (gapped) n-mer features given in D.
% 
%   Example
%     S = 'GTGTCCAGTAGCAAGTTATAAATATCGGACTAGTAATCACTATTAGACA' %use long S
%     M = [0 1; 1 1] 
%     % M = binarySpace(n-1); M = [M true(size(M,1),1)]; %all gapped n-mers
%     D = mappingFoldedKspec(S,M) %calculate background distributions
%     T = {'GTGTCAGTAG','CAAGTTATAA','ATATGGACTA'}
%     F = mappingFunction(T,M,D)
%
%   See also mappingFoldedKspec, binarySpace
if nargin < 4; verbose = false; end
if isempty(M)
    M = binarySpace(length(D(1).model{1})-1); M = [M true(size(M,1),1)];
end
if ischar(T)
    if exist(T,'file')
        FT = fastaread(T);
        T = {};
        for t = 1:size(FT,1)
            T = [T; FT(t).Sequence];
        end
    end
elseif isstruct(T)
    FT = T;
    T = {};
    for t = 1:size(FT,1)
        T = [T; FT(t).Sequence];
    end
end
if iscell(T)
    set_size = max(size(T));
else
    set_size = 1;
    T = {T};
end
if isstruct(D)
    n_models = max(size(D));
else
    n_models = 1;
end
n_models_total = 0;
for i = 1:n_models
    n_models_total = n_models_total + length(D(i).distribution);
end 
F = zeros(set_size,n_models_total);
for t = 1:set_size
    fprintf('Processing sequence %u/%u.. ',t,set_size);
    D_t = mappingFoldedKspec(T{t},M,verbose);
    F_t = [];
    for i = 1:n_models
        F_t = [F_t D_t(i).distribution./D(i).distribution];
    end
    F(t,:) = F_t;
    fprintf('done.\n');
end %t
models = {}; 
for i = 1:n_models
    models = [models D(i).model];
end