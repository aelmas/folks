% ------------------------------------------------------------
n_CVtests = 100;
nFold = 10; 
k = 6; %k-mer
%**************
feature_indices_CV = feature_indices;
%**************
% ------------------------------------------------------------
timeset = fix(clock); testtime = [];
for timelet = timeset 
    testtime = [testtime '-' num2str(timelet)];
end; clear timeset; clear timelet;
% ------------------------------------------------------------
group = [repmat({'negative'},size(negative_sequences,1),1); repmat({'positive'},size(positive_sequences,1),1)];
posclass = 'positive';
% ------------------------------------------------------------
F = [xjn; xip];
% ------------------------------------------------------------
tic;
repeatAUC = [];
gappedkmersFNidx = [];
FPRs = []; TPRs = [];
for repeatCV = 1:n_CVtests
    repeatCV
    trueClasses = [];
    predictedClasses = [];
    AUCs = zeros(nFold,1);
    indices = crossvalind('Kfold',group,nFold)';
    CVidx = [];
    for partition = 1:max(indices)
        testSet = find(indices==partition); %half is positive half is backgr
        trainingSet = find(indices~=partition); %half is posclass half is negclass
        % ------------------------------------------------------------
        myoption = statset('MaxIter',200000);
        svmStruct = svmtrain(F(trainingSet,feature_indices_CV),group(trainingSet),'kernel_function',...
            'linear','options',myoption); %half is posclass half is negclass
        % ------------------------------------------------------------
        predicted_group = [];
        for tS = testSet
            prediction = svmclassify(svmStruct,F(tS,feature_indices_CV));
            predicted_group = [predicted_group; isequal(prediction{1},posclass)];
        end
        [~,~,~,auc_fold] = perfcurve(group(testSet),predicted_group,posclass);
        % ------------------------------------------------------------
        trueClasses = [trueClasses; group(testSet)];
        predictedClasses = [predictedClasses; predicted_group];
        CVidx = [CVidx testSet];
        AUCs(partition) = auc_fold;
    end
    [FPR,TPR,Tscore,AUC] = perfcurve(trueClasses,predictedClasses,posclass);
    meanAUC = mean(AUCs(1:partition));
    repeatAUC = [repeatAUC; AUC];
    for i = 1:size(trueClasses,1)
        if isequal(trueClasses{i},posclass)
            if isequal(predictedClasses(i),0)
                gappedkmersFNidx = [gappedkmersFNidx CVidx(i)];
            end
        end
    end
    figure(12); hold on; plot(FPR,TPR);
    FPRs = [FPRs; FPR]; TPRs = [TPRs; TPR];
end
elapsedtime = toc;
mean(repeatAUC) 
title(sprintf('ROC. AUC_{avg} = %f in %u CV tests with SVM with folded %u-spectrum kernel.',mean(repeatAUC),n_CVtests,k));
ylabel('True positive rate'); xlabel('False positive rate'); box off;
% ------------------------------------------------------------
clear F; save(sprintf('tmpresultsSVMfolded_CV%s.mat',testtime));
% ------------------------------------------------------------
