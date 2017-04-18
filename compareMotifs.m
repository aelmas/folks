function out = compareMotifs(motifA, motifB)
% out = compareMotifs(motifA, motifB) 
% motifA, a query motif
%     example format 1 
%       motifA = ['AACGAA'; 'CACAAA'; 'CACGAA'; 'CCCTAA' ];
%     example format 2
%       motifA = [0.25 0.75 0 0; 0.75 0.25 0 0; 0 1 0 0; 0.25 0 0.5 0.25; 1 0 0 0; 1 0 0 0];
% motifB, database motif, 
%     either in motifA format, or a path to a database file
%       (default) motifB = 'dmel_15_TFBSs_and_JASPAR_CORE_2014_insects.meme'
dbfilename = 'dmel_15_TFBSs_and_JASPAR_CORE_2014_insects.meme';
if ~exist('motifB','var')
    motifB = dbfilename;
elseif isempty(motifB)
    motifB = dbfilename;
end
dbfile = textread(dbfilename,'%s');
alphabet_length = 4;
source_sites_def = 20;
source_Evalue_def = 0;
if ischar(motifA)
    if size(motifA,1)==1
        if isempty(setdiff(upper(motifA),'ACGTN'))
            consensus = upper(motifA);
            motifA = zeros(4,length(consensus));
            for c = 1:length(consensus)
                if isequal(consensus(c),'N')
                    motifA(:,c) = 1/4;
                else
                    motifA(:,c) = 'ACGT'==consensus(c);
                end
            end
            motifA = motifA';
        end
    else
        nucCount = [sum(motifA=='A',1);...
            sum(motifA=='C',1);...
            sum(motifA=='G',1);...
            sum(motifA=='T',1)];
        motifA = (nucCount) ./ repmat(sum(nucCount,1),4,1);
        if ~isempty(find(sum(nucCount,1)==0))
            motifA(:,sum(nucCount,1)==0) = ones(4,length(find(sum(nucCount,1)==0)))/4;
        end
        motifA = motifA';
    end
end
% write motifA
fid = fopen('motifA.txt','w');
fprintf(fid,'MEME version 4\n'); 
fprintf(fid,'ALPHABET= ACGT\n');
fprintf(fid,'MOTIF A\n');
fprintf(fid,'letter-probability matrix: alength= %u w= %u nsites= %u E= %g\n'...
    ,alphabet_length,size(motifA,1),source_sites_def,source_Evalue_def);
for i = 1:size(motifA,1)
    for j = 1:size(motifA,2)
        fprintf(fid,'%g\t',motifA(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
if isnumeric(motifB)
    % write motifB
    fid = fopen('motifB.txt','w');
    fprintf(fid,'MEME version 4\n');
    fprintf(fid,'ALPHABET= ACGT\n');
    fprintf(fid,'MOTIF B\n');
    fprintf(fid,'letter-probability matrix: alength= %u w= %u nsites= %u E= %g\n'...
        ,alphabet_length,size(motifB,1),source_sites_def,source_Evalue_def);
    for i = 1:size(motifB,1)
        for j = 1:size(motifB,2)
            fprintf(fid,'%g\t',motifB(i,j));
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    % perform tomtom
%     [~,result] = 
    system(['./tomtom -oc tomtom_results motifA.txt motifB.txt']);
elseif ischar(motifB)
%     [~,result] = 
    system(['./tomtom -oc tomtom_results motifA.txt ' motifB]);
end
output = importdata('tomtom_results/tomtom.txt');
out.Header = textscan(output{1},'%s'); out.Header = out.Header{1};
for hit = 1:size(output,1)-1
    out.Hit{hit} = textscan(output{hit+1},'%s'); 
    for i = 1:length(dbfile)
        if isequal(dbfile{i},out.Hit{hit}{1}{2})
            out.Hit{hit}{1} = [out.Hit{hit}{1}; dbfile{i+1}];
            break
        end
    end
end

