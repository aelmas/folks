load Stringham_2013_CRMSwithFlanking_filenames.mat
dataset = 'trueCRMs';
filenames = filenames0; CRMxmap = 1:size(filenames,1);

fid = fopen(sprintf('tmp_results_%s.xls',dataset),'w');          %*****
fprintf(fid,'Table 1. Top enriched features (weight cutoff %g) in each CRM. SVM is trained by the entire set (%u CRMs, %u background)\n',r_cutoff,size(positive_sequences,1),size(negative_sequences,1));
fprintf(fid,'Final number of features after feature elimination = %u\n',length(feature_indices));
%    
fprintf(fid,'GROUP 1: The group of %u CRMs that are detected by gapped k-mers but missed by contiguous k-mers\n',length(Group1_CRMs));
fprintf(fid,'GROUP 2: The group of %u CRMs that are missed by both gapped k-mers and contiguous k-mers\n',length(Group2_CRMs));
fprintf(fid,'GROUP 3: The group of %u CRMs that are detected by contiguous k-mers but missed by gapped k-mers\n',length(Group3_CRMs));
fprintf(fid,'GROUP 4: The group of %u CRMs that are detected by both contiguous k-mers and gapped k-mers\n',length(Group4_CRMs));
fprintf(fid,'\t');
for i = 1:length(Group1_CRMs); fprintf(fid,'GROUP 1\t\t'); end
for i = 1:length(Group2_CRMs); fprintf(fid,'GROUP 2\t\t'); end
for i = 1:length(Group3_CRMs); fprintf(fid,'GROUP 3\t\t'); end
for i = 1:length(Group4_CRMs); fprintf(fid,'GROUP 4\t\t'); end
%
fprintf(fid,'\n'); fprintf(fid,'\t');
for i = 1:length(CRM_indices); fprintf(fid,'CRM ID: %u\t\t',CRMxmap(CRM_indices(i))); end
fprintf(fid,'\n');
for c = 1:length(CRM_indices)
    fprintf(fid,'\t%s\tWeight',filenames{CRM_indices(c)});
end; fprintf(fid,'\n');
for r = 1:maxFeaturesInCRM
    fprintf(fid,'%u\t',r);
    for c = 1:length(CRM_indices)
        if r <= size(enriched_features{c},1)
            fprintf(fid,'%s\t',enriched_features{c}(r,:));
            fprintf(fid,'%12.5f\t',enriched_features_weights{c}(r));
        else
            fprintf(fid,'\t\t');
        end
    end
    fprintf(fid,'\n');
end
%
fprintf(fid,'\n\nTable 2. Corresponding gapped k-mers\n');
fprintf(fid,'\t');
for i = 1:length(CRM_indices); fprintf(fid,'CRM ID: %u\t\t',CRMxmap(CRM_indices(i))); end
fprintf(fid,'\n');
for c = 1:length(CRM_indices)
    fprintf(fid,'\t%s\tWeight',filenames{CRM_indices(c)});
end; fprintf(fid,'\n');
for r = 1:maxFeaturesInCRM
    fprintf(fid,'%u\t',r);
    for c = 1:length(CRM_indices)
        if r <= size(enriched_gapped_kmers{c},1)
            fprintf(fid,'%s\t',enriched_gapped_kmers{c}(r,:));
            fprintf(fid,'%12.5f\t',enriched_features_weights{c}(r));
        else
            fprintf(fid,'\t\t');
        end
    end
    fprintf(fid,'\n');
end
%
all_different_gapped_kmers = [];
all_different_gapped_kmers_counts = [];
different_gapped_kmers = cell(1,length(CRM_indices));
different_gapped_kmers_counts = cell(1,length(CRM_indices));
different_gapped_kmers_weights = cell(1,length(CRM_indices));
for c = 1:length(CRM_indices)
    rank_data_indices = [];
    for r = 1:size(enriched_gapped_kmers{c},1)
        is_counted = false;
        for g = 1:size(different_gapped_kmers{c},1)
            if isequal(enriched_gapped_kmers{c}(r,:),different_gapped_kmers{c}(g,:))
                is_counted = true;
                break
            end
        end
        if ~is_counted
            count = 0;
            total_weight = 0;
            for r2 = 1:size(enriched_gapped_kmers{c},1)
                if isequal(enriched_gapped_kmers{c}(r,:),enriched_gapped_kmers{c}(r2,:))
                    count = count+1;
                    total_weight = total_weight + enriched_features_weights{c}(r2);
                end
            end
            all_different_gapped_kmers = [all_different_gapped_kmers; enriched_gapped_kmers{c}(r,:)]; 
            all_different_gapped_kmers_counts = [all_different_gapped_kmers_counts; count]; 
            different_gapped_kmers{c} = [different_gapped_kmers{c}; enriched_gapped_kmers{c}(r,:)];
            different_gapped_kmers_counts{c} = [different_gapped_kmers_counts{c}; count];
            different_gapped_kmers_weights{c} = [different_gapped_kmers_weights{c}; total_weight];
            query = enriched_gapped_kmers{c}(r,:);
            while isequal(lower(query(1)),'n'); query = query(2:end); end
            while isequal(lower(query(end)),'n'); query = query(1:end-1); end
        end
    end %r
end %c
%
unique_gapped_kmers = [];
unique_gapped_kmers_counts = [];
unique_gapped_kmers_appears_in_CRMS = [];
for r = 1:size(all_different_gapped_kmers,1)
    is_counted = false;
    for g = 1:size(unique_gapped_kmers,1)
        if isequal(all_different_gapped_kmers(r,:),unique_gapped_kmers(g,:))
            unique_gapped_kmers_counts(g) = unique_gapped_kmers_counts(g) + all_different_gapped_kmers_counts(r);
            unique_gapped_kmers_appears_in_CRMS(g) = unique_gapped_kmers_appears_in_CRMS(g) + 1;    
            is_counted = true;
            break
        end
    end
    if ~is_counted
        unique_gapped_kmers = [unique_gapped_kmers; all_different_gapped_kmers(r,:)];
        unique_gapped_kmers_counts = [unique_gapped_kmers_counts; all_different_gapped_kmers_counts(r)];
        unique_gapped_kmers_appears_in_CRMS = [unique_gapped_kmers_appears_in_CRMS; 1];
    end
end
%
maxNumDifGkmers = 1;
for c = 1:length(CRM_indices)
    if maxNumDifGkmers < size(different_gapped_kmers_counts{c},1)
        maxNumDifGkmers = size(different_gapped_kmers_counts{c},1);
    end
end
%
fprintf(fid,'\n\nTable 3. Counts and cumulative weights of different gapped k-mers. (Rank-ordered by weights as in the Tables 1-2)\n');
fprintf(fid,'\t');
for i = 1:length(CRM_indices); fprintf(fid,'CRM ID: %u\t\t',CRMxmap(CRM_indices(i))); end
fprintf(fid,'\n');
for c = 1:length(CRM_indices)
    fprintf(fid,'\t%s: Count\tCumulative Weight',filenames{CRM_indices(c)});
end; fprintf(fid,'\n');
for r = 1:maxNumDifGkmers
    fprintf(fid,'%u\t',r);
    for c = 1:length(CRM_indices)
        if r <= size(different_gapped_kmers_counts{c},1)
            fprintf(fid,'%s: %u\t',different_gapped_kmers{c}(r,:),different_gapped_kmers_counts{c}(r));
            fprintf(fid,'%12.5f\t',different_gapped_kmers_weights{c}(r));
        else
            fprintf(fid,'\t\t');
        end
    end
    fprintf(fid,'\n');
end
%
fprintf(fid,'\n\nTable 4. Total number of occurrences of each different gapped k-mer (TotalCount) across the CRMs, and its relative occurrences, i.e., Count/(TotalCount/#CRMs the kmer found in).\t \n');
fprintf(fid,'\t');
for i = 1:length(CRM_indices); fprintf(fid,'CRM ID: %u\t\t',CRMxmap(CRM_indices(i))); end
fprintf(fid,'\n');
for c = 1:length(CRM_indices)
    fprintf(fid,'\t%s: TotalCount in #CRMs\tCount/(TotalCount/#CRMs)',filenames{CRM_indices(c)});
end; fprintf(fid,'\n');
different_gapped_kmers_CountRatios = cell(1,length(CRM_indices));
for r = 1:maxNumDifGkmers+1
    if r <= maxNumDifGkmers
        fprintf(fid,'%u\t',r);
    else
        fprintf(fid,'\t');
    end
    for c = 1:length(CRM_indices)
        if r <= size(different_gapped_kmers_counts{c},1)
            total_count = 0;
            totol_appearance = 0;
            for i = 1:size(unique_gapped_kmers,1)
                if isequal(different_gapped_kmers{c}(r,:),unique_gapped_kmers(i,:))
                    total_count = unique_gapped_kmers_counts(i);
                    total_appearance = unique_gapped_kmers_appears_in_CRMS(i);
                    break;
                end
            end
            count_ratio = different_gapped_kmers_counts{c}(r)/(total_count/total_appearance);
            fprintf(fid,'%s: %u in %u\t',different_gapped_kmers{c}(r,:),total_count,total_appearance);
            fprintf(fid,'%12.4f\t',count_ratio);
            different_gapped_kmers_CountRatios{c} = [different_gapped_kmers_CountRatios{c}; count_ratio];
        elseif r == maxNumDifGkmers+1
            if isempty(different_gapped_kmers_counts{c});
                fprintf(fid,'\t\t');
            else
                fprintf(fid,'average count ratio\t%12.4f\t',mean(different_gapped_kmers_CountRatios{c}));
            end
        else
            fprintf(fid,'\t\t');
        end
    end
    fprintf(fid,'\n');
end

