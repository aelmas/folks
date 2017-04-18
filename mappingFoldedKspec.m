function D = mappingFoldedKspec(S,M,verbose)
% mappingFoldedKspec Computes the distribution of given sequence model(s)
%   D = mappingFoldedKspec(S,M,verbose) calculates the distribution of each
%   entries (rows) of the sequence models M in the sequence S.
%
%   S, a character string representing a DNA sequence, or is the path to
%     the fasta file containing a DNA sequence.
%
%   M, a matrix where each row represents a gapped n-mer model.
%   Example
%     [1 1 1] represents all 3-mers ['AAA','AAC',...,'ACA','ACC',...,'TTT']
%     [1 0 1] represents the gapped 3-mers where the 2nd nucleotide
%       is ignored ['AgA','AgC',...,'CgA','CgC',...,'TgT']
%     Note: [1 0 1] have the same number of elements with 2-mers.
%
%   verbose, logical value for displaying the progress. default is false.
%
%   D, struct array. Each row contains the corresponding (gapped) n-mer
%     distributions.
%
%   Example
%     S = 'AAATGAATTGC'
%     M = [0 1; 1 1]
%     % M = binarySpace(n-1); M = [M true(size(M,1),1)]; %all gapped n-mers
%     D = mappingFoldedKspec(S,M)
%       D(1) returns the mononucleotide distributions due to [0 1] of M
%           distribution    model
%                   4/10     'gA'
%                   1/10     'gC'
%                   2/10     'gG'
%                   3/10     'gT'
%       D(2) returns the dinonucleotide distributions due to [1 1] of M
%           distribution    model
%                    3/9     'AA'
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
%   See also mappingKspec, binarySpace
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
    end
end
clear FS;
[sM1,sM2] = size(M);
alphabet = 'ACGT'; l = length(alphabet);
% compute contiguous n-mers distribution
if isempty(find(sum(M,2)==size(M,2)))
    D_nmer = mappingKspec(S,ones(1,sM2));
else
    D_nmer = mappingKspec(S,ones(1,sM2),verbose);
end
% compute gapped n-mers distribution
% D = [];
for i = 1:sM1
    ms = sum(M(i,:));
    if verbose
        fprintf('model %u/%u : %s .. \n',i,sM1,mat2str(double(M(i,:))));
    end
    if ms == sM2 %the longest contiguous k-mer
        D(i).distribution = D_nmer.distribution;
        D(i).model = D_nmer.model;
        if verbose; fprintf(' done.\n'); end
        continue
    end        
    D(i).distribution = zeros(1,l^ms);
    for j = 1:l^ms
        mdl = repmat('g',1,sM2);
        lbase = [];
        for st = dec2base(j-1,l,ms)
            lbase = [lbase str2double(st)];
        end
        mdl(M(i,:)) = alphabet(lbase+1);
        D(i).model(j) = {mdl};
        if verbose
            fprintf('%u/%u: %s .. ',j,l^ms,mdl);
        end
        %
        myindex = [];
        for k = 1:l^(sM2-ms)
            lgap = [];
            for st = dec2base(k-1,l,sM2-ms)
                lgap = [lgap str2double(st)];
            end
            binbase = zeros(1,sM2);
            binbase(M(i,:)) = lbase;
            binbase(~M(i,:)) = lgap;
            myindex = [myindex 1+4.^(sM2-1:-1:0)*(binbase)'];
        end
        D(i).distribution(j) = sum(D_nmer.distribution(myindex));
        %
        if verbose; fprintf(' done.\n'); end
    end
    if verbose; fprintf('done.\n'); end
end

end %function
%