function space = binarySpace(n)
% binarySpace Output the sample space consisting of all binary vectors of 
% a given length
%   binarySpace(n) returns the sample space of binary vectors (rows)
%   of length-n in the order of increasing decimal values
%
%   Example
%       binarySpace(2) returns [0 0; 0 1; 1 0; 1 1] in logical matrix
%
%   See also dec2bin, bin2dec, mappingFoldedKspec
space = false(2^n,n);
for i = 1:2^n
    tmp = dec2bin(i-1,n);
    for j = 1:n
        if isequal(tmp(j),'1')
            space(i,j) = true;
        end
    end
end