function space = binarySpace(k)
% binarySpace Output the sample space consisting of all binary vectors of 
% a given length
%   binarySpace(k) returns the sample space of binary vectors (rows)
%   of length-k in the order of increasing decimal values
%
%   Example
%       binarySpace(2) returns [0 0; 0 1; 1 0; 1 1] in logical matrix
%
%   See also dec2bin, bin2dec, mappingFoldedKspec
space = false(2^k,k);
for i = 1:2^k
    tmp = dec2bin(i-1,k);
    for j = 1:k
        if isequal(tmp(j),'1')
            space(i,j) = true;
        end
    end
end