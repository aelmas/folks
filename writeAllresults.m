load Stringham_2013_CRMSwithFlanking_filenames.mat
dataset = 'trueCRMs';
filenames = filenames0; CRMxmap = 1:size(filenames,1);

fid = fopen(sprintf('tmp_results_%s.xls',dataset),'w');          %*****
fprintf(fid,'Table 1. Top enriched features (weight cutoff %g) in each CRM. SVM is trained by the entire set (%u CRMs, %u background)\n',r_cutoff,size(positive_sequences,1),size(negative_sequences,1));
fprintf(fid,'Final number of features after feature elimination = %u\n',length(feature_indices));
%    
fprintf(fid,'GROUP 1: The group of %u CRMs that are detected by gapped n-mers but missed by contiguous n-mers\n',length(Group1_CRMs));
fprintf(fid,'GROUP 2: The group of %u CRMs that are missed by both gapped n-mers and contiguous n-mers\n',length(Group2_CRMs));
fprintf(fid,'GROUP 3: The group of %u CRMs that are detected by contiguous n-mers but missed by gapped n-mers\n',length(Group3_CRMs));
fprintf(fid,'GROUP 4: The group of %u CRMs that are detected by both contiguous n-mers and gapped n-mers\n',length(Group4_CRMs));
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
fprintf(fid,'\n\nTable 2. Corresponding gapped n-mers\n');
fprintf(fid,'\t');
for i = 1:length(CRM_indices); fprintf(fid,'CRM ID: %u\t\t',CRMxmap(CRM_indices(i))); end
fprintf(fid,'\n');
for c = 1:length(CRM_indices)
    fprintf(fid,'\t%s\tWeight',filenames{CRM_indices(c)});
end; fprintf(fid,'\n');
for r = 1:maxFeaturesInCRM
    fprintf(fid,'%u\t',r);
    for c = 1:length(CRM_indices)
        if r <= size(enriched_gapped_nmers{c},1)
            fprintf(fid,'%s\t',enriched_gapped_nmers{c}(r,:));
            fprintf(fid,'%12.5f\t',enriched_features_weights{c}(r));
        else
            fprintf(fid,'\t\t');
        end
    end
    fprintf(fid,'\n');
end
%
all_different_gapped_nmers = [];
all_different_gapped_nmers_counts = [];
different_gapped_nmers = cell(1,length(CRM_indices));
different_gapped_nmers_counts = cell(1,length(CRM_indices));
different_gapped_nmers_weights = cell(1,length(CRM_indices));
for c = 1:length(CRM_indices)
    rank_data_indices = [];
    for r = 1:size(enriched_gapped_nmers{c},1)
        is_counted = false;
        for g = 1:size(different_gapped_nmers{c},1)
            if isequal(enriched_gapped_nmers{c}(r,:),different_gapped_nmers{c}(g,:))
                is_counted = true;
                break
            end
        end
        if ~is_counted
            count = 0;
            total_weight = 0;
            for r2 = 1:size(enriched_gapped_nmers{c},1)
                if isequal(enriched_gapped_nmers{c}(r,:),enriched_gapped_nmers{c}(r2,:))
                    count = count+1;
                    total_weight = total_weight + enriched_features_weights{c}(r2);
                end
            end
            all_different_gapped_nmers = [all_different_gapped_nmers; enriched_gapped_nmers{c}(r,:)]; 
            all_different_gapped_nmers_counts = [all_different_gapped_nmers_counts; count]; 
            different_gapped_nmers{c} = [different_gapped_nmers{c}; enriched_gapped_nmers{c}(r,:)];
            different_gapped_nmers_counts{c} = [different_gapped_nmers_counts{c}; count];
            different_gapped_nmers_weights{c} = [different_gapped_nmers_weights{c}; total_weight];
            query = enriched_gapped_nmers{c}(r,:);
            while isequal(lower(query(1)),'g'); query = query(2:end); end
            while isequal(lower(query(end)),'g'); query = query(1:end-1); end
        end
    end %r
end %c
%
unique_gapped_nmers = [];
unique_gapped_nmers_counts = [];
unique_gapped_nmers_appears_in_CRMS = [];
for r = 1:size(all_different_gapped_nmers,1)
    is_counted = false;
    for g = 1:size(unique_gapped_nmers,1)
        if isequal(all_different_gapped_nmers(r,:),unique_gapped_nmers(g,:))
            unique_gapped_nmers_counts(g) = unique_gapped_nmers_counts(g) + all_different_gapped_nmers_counts(r);
            unique_gapped_nmers_appears_in_CRMS(g) = unique_gapped_nmers_appears_in_CRMS(g) + 1;    
            is_counted = true;
            break
        end
    end
    if ~is_counted
        unique_gapped_nmers = [unique_gapped_nmers; all_different_gapped_nmers(r,:)];
        unique_gapped_nmers_counts = [unique_gapped_nmers_counts; all_different_gapped_nmers_counts(r)];
        unique_gapped_nmers_appears_in_CRMS = [unique_gapped_nmers_appears_in_CRMS; 1];
    end
end
%
maxNumDifGnmers = 1;
for c = 1:length(CRM_indices)
    if maxNumDifGnmers < size(different_gapped_nmers_counts{c},1)
        maxNumDifGnmers = size(different_gapped_nmers_counts{c},1);
    end
end
%
fprintf(fid,'\n\nTable 3. Counts and cumulative weights of different gapped n-mers. (Rank-ordered by weights as in the Tables 1-2)\n');
fprintf(fid,'\t');
for i = 1:length(CRM_indices); fprintf(fid,'CRM ID: %u\t\t',CRMxmap(CRM_indices(i))); end
fprintf(fid,'\n');
for c = 1:length(CRM_indices)
    fprintf(fid,'\t%s: Count\tCumulative Weight',filenames{CRM_indices(c)});
end; fprintf(fid,'\n');
for r = 1:maxNumDifGnmers
    fprintf(fid,'%u\t',r);
    for c = 1:length(CRM_indices)
        if r <= size(different_gapped_nmers_counts{c},1)
            fprintf(fid,'%s: %u\t',different_gapped_nmers{c}(r,:),different_gapped_nmers_counts{c}(r));
            fprintf(fid,'%12.5f\t',different_gapped_nmers_weights{c}(r));
        else
            fprintf(fid,'\t\t');
        end
    end
    fprintf(fid,'\n');
end
%
fprintf(fid,'\n\nTable 4. Total number of occurrences of each different gapped n-mer (TotalCount) across the CRMs, and its relative occurrences, i.e., Count/(TotalCount/#CRMs the nmer found in).\t \n');
fprintf(fid,'\t');
for i = 1:length(CRM_indices); fprintf(fid,'CRM ID: %u\t\t',CRMxmap(CRM_indices(i))); end
fprintf(fid,'\n');
for c = 1:length(CRM_indices)
    fprintf(fid,'\t%s: TotalCount in #CRMs\tCount/(TotalCount/#CRMs)',filenames{CRM_indices(c)});
end; fprintf(fid,'\n');
different_gapped_nmers_CountRatios = cell(1,length(CRM_indices));
for r = 1:maxNumDifGnmers+1
    if r <= maxNumDifGnmers
        fprintf(fid,'%u\t',r);
    else
        fprintf(fid,'\t');
    end
    for c = 1:length(CRM_indices)
        if r <= size(different_gapped_nmers_counts{c},1)
            total_count = 0;
            totol_appearance = 0;
            for i = 1:size(unique_gapped_nmers,1)
                if isequal(different_gapped_nmers{c}(r,:),unique_gapped_nmers(i,:))
                    total_count = unique_gapped_nmers_counts(i);
                    total_appearance = unique_gapped_nmers_appears_in_CRMS(i);
                    break;
                end
            end
            count_ratio = different_gapped_nmers_counts{c}(r)/(total_count/total_appearance);
            fprintf(fid,'%s: %u in %u\t',different_gapped_nmers{c}(r,:),total_count,total_appearance);
            fprintf(fid,'%12.4f\t',count_ratio);
            different_gapped_nmers_CountRatios{c} = [different_gapped_nmers_CountRatios{c}; count_ratio];
        elseif r == maxNumDifGnmers+1
            if isempty(different_gapped_nmers_counts{c});
                fprintf(fid,'\t\t');
            else
                fprintf(fid,'average count ratio\t%12.4f\t',mean(different_gapped_nmers_CountRatios{c}));
            end
        else
            fprintf(fid,'\t\t');
        end
    end
    fprintf(fid,'\n');
end

