% =========================================================================
% SVM with Feature Elimination (basic)
%   variable names follow the same notations as in the paper
% =========================================================================
clear, clc
r_cutoff = 0.005; k = 6;
% -------------------------------------------------------------------------
% L. load sequence data 
positive_sequences = fastaread('true_positive_set_Stringham13flanking.fasta');
% -------------------------------------------------------------------------
% 0. generate positive training data (x_i^p)
M = binarySpace(k-1); M = [M true(size(M,1),1)]; %all gapped k-mers
%     load dmel_genome_ucsc_2006_dm3
%     backgroundseq = []; 
%     for i = 10 %size(dmel_genome_ucsc_2006_dm3,1)
%         backgroundseq = [backgroundseq upper(dmel_genome_ucsc_2006_dm3(i).Sequence)]; 
%     end; clear dmel_genome_ucsc_2006_dm3
%     D = mappingFoldedKspec(backgroundseq,M); %generate background frequencies of gapped k-mers 
load D_dmel_genomeAll_gapped6mers %load background frequencies of gapped 6-mers in Drosophila Melanogaster genome
xip = mappingFunction(upper(positive_sequences),M,D);
% -------------------------------------------------------------------------
% 1. generate negative training data (x_j^n)
negative_sequences = {};
NSscramble = 10; 
for i = 1:size(positive_sequences,1)
    tmpseq = positive_sequences(i).Sequence;
    for j = 1:NSscramble
        negative_sequences = [negative_sequences; tmpseq( randperm(length(tmpseq)) ) ];
    end
end
[xjn,features] = mappingFunction(upper(negative_sequences),M,D); 
% -------------------------------------------------------------------------
% 2. generate false positive training data (x_j^f)
false_positive_sequences = {};
FPSscramble = 10; 
for i = 1:size(positive_sequences,1)
    tmpseq = positive_sequences(i).Sequence;
    for j = 1:FPSscramble
        false_positive_sequences = [false_positive_sequences; tmpseq( randperm(length(tmpseq)) ) ];
    end
end
xjf = mappingFunction(upper(false_positive_sequences),M,D); 
% =========================================================================
% FEATURE ELIMINATION
% =========================================================================
% 3. train SVM with (-,+) data
svm_settings = statset('MaxIter',3e5);
svmStruct = svmtrain([xjn;xip],[repmat('-',size(xjn,1),1);repmat('+',size(xip,1),1)],...
    'kernel_function','linear','options',svm_settings);
% -------------------------------------------------------------------------
% 3f. train SVM with (-,+f) data
svmStructf = svmtrain([xjn;xjf],[repmat('-',size(xjn,1),1);repmat('+',size(xjf,1),1)],...
    'kernel_function','linear','options',svm_settings);
% -------------------------------------------------------------------------
% 4. calculate decision vector w for step 3
w = zeros(1,length(svmStruct.SupportVectors(1,:)));
for i = 1:length(svmStruct.SupportVectorIndices)
    group_sign = -1; %+1 first group (background), -1 second group (sequence)
    w = w + svmStruct.Alpha(i)*group_sign * svmStruct.SupportVectors(i,:);
end
% -------------------------------------------------------------------------
% 4f. calculate decision vector wf for step 3f
wf = zeros(1,length(svmStructf.SupportVectors(1,:)));
for i = 1:length(svmStructf.SupportVectorIndices)
    group_sign = -1; %+1 first group (background), -1 second group (sequence)
    wf = wf + svmStructf.Alpha(i)*group_sign * svmStructf.SupportVectors(i,:);
end
% -------------------------------------------------------------------------
% 5. calculate enrichments for each input sequence
ri = xip .* repmat(w,size(xip,1),1);
% -------------------------------------------------------------------------
% 5f. calculate false enrichments for each input sequence
rjf = xjf .* repmat(wf,size(xjf,1),1);
% -------------------------------------------------------------------------
% 6. filter out r_i <= r_cutoff
ri(ri<=r_cutoff) = 0;
% -------------------------------------------------------------------------
% 6f. filter out r_j^f <= r_cutoff
rjf(rjf<=r_cutoff) = 0;
% -------------------------------------------------------------------------
% 7. FE: eliminate r_j^f features from the corresponding r_i
feature_indices = [];
for i = 1:size(ri,1)
    features_i = find(ri(i,:)>r_cutoff);
    for j = (i-1)*10+(1:10)        
        features_i = setdiff(features_i,find(rjf(j,:)>r_cutoff));                
    end
    feature_indices = [feature_indices features_i];
end
% -------------------------------------------------------------------------
% 8. constrain the feature set to the union of remaining features
feature_indices = union(feature_indices,feature_indices);
% =========================================================================
% FINAL PREDICTION
% =========================================================================
% 9. train SVM with (-,+) data
svmStructh = svmtrain([xjn(:,feature_indices);xip(:,feature_indices)],[repmat('-',size(xjn,1),1);repmat('+',size(xip,1),1)],...
    'kernel_function','linear','options',svm_settings);
% -------------------------------------------------------------------------
% 10. calculate decision vector wh for step 9
wh = zeros(1,length(svmStructh.SupportVectors(1,:)));
for i = 1:length(svmStructh.SupportVectorIndices)
    group_sign = -1; %+1 first group (background), -1 second group (sequence)
    wh = wh + svmStructh.Alpha(i)*group_sign * svmStructh.SupportVectors(i,:);
end
% -------------------------------------------------------------------------
% 11. calculate final enrichments for each input sequence
ri_final = xip(:,feature_indices) .* repmat(wh,size(xip,1),1);
% -------------------------------------------------------------------------
% 12. filter out r_i_final <= r_cutoff
ri_final(ri_final<=r_cutoff) = 0;
% -------------------------------------------------------------------------
% RESULTS: features_final, ri_final
features_final = {};
for i = 1:size(ri_final,1)
    features_final{i} = features(ri_final(i,:)>r_cutoff)';
end


%% WRITE ENRICHMENT TABLES
pipelineSVMCVregular
pipelineSVMCVfolded

n_CVtests = 100;
minFNrate = .1;
kmersFNidx0 = kmersFNidx;
kmersFNidx = [];
for i = uniondata(kmersFNidx0)
    if length(find(kmersFNidx0==i)) > minFNrate * n_CVtests
        kmersFNidx = [kmersFNidx kmersFNidx0(kmersFNidx0==i)];
    end
end
gappedkmersFNidx0 = gappedkmersFNidx;
gappedkmersFNidx = [];
for i = uniondata(gappedkmersFNidx0)
    if length(find(gappedkmersFNidx0==i)) > minFNrate * n_CVtests
        gappedkmersFNidx = [gappedkmersFNidx gappedkmersFNidx0(gappedkmersFNidx0==i)];
    end
end

controlSize = size(negative_sequences,1);
signalSize = size(positive_sequences,1);
Group1_CRMs = setdiff(kmersFNidx,gappedkmersFNidx)-controlSize
Group2_CRMs = intersect(kmersFNidx,gappedkmersFNidx)-controlSize
Group3_CRMs = setdiff(gappedkmersFNidx,kmersFNidx)-controlSize
Group4_CRMs = setdiff(1:signalSize,[Group1_CRMs Group2_CRMs Group3_CRMs]);
CRM_indices = [Group1_CRMs Group2_CRMs Group3_CRMs Group4_CRMs];

xS = 1;
kmerfncount = zeros(1,length(CRM_indices));
gappedkmerfncount = zeros(1,length(CRM_indices));
for j = 1:length(CRM_indices)
    range = (j-1)*xS+1:j*xS;
    for i = 1:length(kmersFNidx)
        if range(1) <= kmersFNidx(i)-controlSize && kmersFNidx(i)-controlSize <= range(end)
            kmerfncount(j) = kmerfncount(j) + 1;
        end
    end
    for i = 1:length(gappedkmersFNidx)
        if range(1) <= gappedkmersFNidx(i)-controlSize && gappedkmersFNidx(i)-controlSize <= range(end)
            gappedkmerfncount(j) = gappedkmerfncount(j) + 1;
        end
    end
end
kmerfncount = kmerfncount/xS;
gappedkmerfncount = gappedkmerfncount/xS;

    figure; bar([kmerfncount' gappedkmerfncount'])
    box off; axis([1 length(kmerfncount)+1 0 n_CVtests]); grid on;
    xlabel('CRM ID'); ylabel('# FN-predicted tests');
    legend('6-spectrum kernel (FE)','Folded 6-spectrum kernel (FE)',...
        'Location','northoutside','Orientation','horizontal'); legend('boxoff')
    fig = gcf; 
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 2];
    fig.PaperPositionMode = 'manual';
    print(sprintf('tmp_plot_FNnumbers'),'-depsc','-tiff','-r300','-loose')
    
full_feature_set = 1:size(features,2);    
enriched_features = cell(1,length(CRM_indices));
enriched_features_weights = cell(1,length(CRM_indices));
enriched_gapped_kmers = cell(1,length(CRM_indices));
maxFeaturesInCRM = 0;
for c = 1:length(CRM_indices);
    wfx = ri(CRM_indices(c),:);
    [wfx_sorted,sidx] = sort(wfx,'descend');
    feature_indices_sorted = full_feature_set(sidx);
    %
    idxFeaturesThresh_sorted = feature_indices_sorted(wfx_sorted>r_cutoff);   
    %
    idxFeaturesThresh_sorted = intersect(idxFeaturesThresh_sorted, features_i_h{CRM_indices(c)},'stable');
    % ---------------------------------------------------------------
    enriched_features{c} = cell2mat(features(idxFeaturesThresh_sorted)');
    enriched_features_weights{c} = wfx(idxFeaturesThresh_sorted);
    % ---------------------------------------------------------------
    if length(enriched_features_weights{c}) > maxFeaturesInCRM
        maxFeaturesInCRM = length(enriched_features_weights{c});
    end
    for mdl = enriched_features{c}'
        gapped_kmer = repmat('a',1,length(mdl));
        gapped_kmer(upper(mdl)=='N') = 'N';
        enriched_gapped_kmers{c} = [enriched_gapped_kmers{c}; gapped_kmer];
    end
end
all_enriched_features = [];
for i = 1:size(enriched_features,2)
    all_enriched_features = [all_enriched_features; enriched_features{i}];
end
all_unique_enriched_features = uniondata(all_enriched_features')';
%
mmask = zeros(size(features));
for i = 1:size(all_unique_enriched_features,1)
    mmask = mmask + strcmp(all_unique_enriched_features(i,:),features);
end

display('writing..');
writeAllresults

