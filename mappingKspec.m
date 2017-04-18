function D = mappingKspec(S,M,verbose)
% mappingKspec Computes the distribution (frequency) of the given 
%   contagious sequence model(s)
%   D = mappingKspec(S,M) calculates the distribution of the sequence 
%   model M in the sequence S.
%
%   S, a character string representing a DNA sequence, or is the path to 
%     the fasta file containing a DNA sequence.
%
%   M, is a sequence of 1s representing a contiguous n-mer.
%   Example 
%     [1 1 1]   represents all 3-mers ['AAA','AAC',...,'ACA','ACC',...,'TTT']
%     [1 1]     represents all dimers
%
%   verbose, logical value for displaying the progress. default is false.
%
%   D, struct array. [model] Each row contains the corresponding sequence 
%     features, e.g., n-mer distributions. [distribution] Contains the
%     corresponding distribution of the model.
%
%   Example
%     S = 'AATGAATTGC'
%     M = [1 1] 
%     D = mappingKspec(S,M)
%       D returns the dinonucleotide distributions due to [1 1] of M 
%           distribution    model
%                    2/9     'AA'
%                      0     'AC'
%                      0     'AG'
%                    2/9     'AT'
%                      0     'CA'
%                      0     'CC'
%                      0     'CG'
%                      0     'CT'
%                    1/9     'GA'
%                    1/9     'GC'
%                      0     'GG'
%                      0     'GT'
%                      0     'TA'
%                      0     'TC'
%                    2/9     'TG'
%                    1/9     'TT'
%
%   See also binarySpace
if nargin < 3
    verbose = false;
end
if ~isequal(class(M),'logical')
    M = logical(M);
end
if ischar(S)
    if exist(S,'file')
        FS = fastaread(S);
        S = [];
        for s = 1:size(FS,1)
            S = [S FS(s).Sequence];
        end
        clear FS;
    end
end
lS = length(S); 
lSnuc = lS-length(find(S=='N'))-length(find(S=='n')); 
sM2 = size(M,2);
alphabet = 'ACGT'; l = length(alphabet);
for i = 1:size(M,1)
    ms = sum(M(i,:));
    if verbose
        fprintf('Computing %u-mer model : %s .. \n',size(M,2),mat2str(double(M(i,:))));
    end
    D(i).distribution = zeros(1,l^ms);
    %
    for n = 1:lS-sM2+1
        base = S(n:n+sM2-1);
        binbase = zeros(1,sM2); %A
        binbase(base=='C') = 1; %C
        binbase(base=='G') = 2; %G
        binbase(base=='T') = 3; %T
        jstar = 1+4.^(sM2-1:-1:0)*(binbase)';
        D(i).distribution(jstar) = D(i).distribution(jstar) + 1;
    end
    D(i).distribution = D(i).distribution/(lSnuc-sM2+1);
    %
    for j = 1:l^ms
        mdl = repmat('g',1,sM2); 
        lbase = dec2base(j-1,l,ms);
        msites = []; 
        for m = lbase
            msites = [msites alphabet(str2double(m)+1)];
        end
        mdl(M(i,:)) = msites;
        if verbose
            fprintf('%u/%u: %s .. ',j,l^ms,mdl);
        end
        D(i).model(j) = {mdl};
        if verbose; fprintf(' done.\n'); end
    end
    if verbose; fprintf('done.\n'); end
end
