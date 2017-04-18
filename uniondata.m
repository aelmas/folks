function [U,C,I] = uniondata(X)
% [unified_columns,counts,indices_in_X] = uniondata(X)
ncols = size(X,2);
if ncols == 0
    U = X; %Unified columns
    C = 1; %Counts
    I = 1;
    return
else
    U = X; 
    U(:,1) = X(:,1);
    C = zeros(1,ncols); 
    C(1) = 1;
    I = zeros(1,ncols);
    I(1) = 1;
    nUnique = 1;
    for i = 2:ncols
        unique = true;
        for j = 1:nUnique
            if isequal(X(:,i),U(:,j))
                C(j) = C(j)+1;
                unique = false;
                break
            end
        end
        if unique
            nUnique = nUnique+1;
            U(:,nUnique) = X(:,i);
            C(nUnique) = 1;
            I(nUnique) = i;
        end
    end
    U = U(:,1:nUnique);
    C = C(1:nUnique);
    I = I(1:nUnique);
end
end