%%
extracted_labels = cellfun(@(x) regexprep(x, '(\w+\s+\w+).*', '$1'), categories_samples, 'UniformOutput', false);
%%
% Inside the loop (after training)
YPred = cell(1, kfoldsize);
scores = cell(1, kfoldsize);
cm_s = cell(1, kfoldsize);
cm_sum = zeros(length(extracted_labels), length(extracted_labels));
YValidation = cell(1, kfoldsize);


%%
for n = 1:kfoldsize
    YValidation{n} = validation_cat{n};
    [YPred{n}, scores{n}] = classify(vnet{n}, validation_dados{n});
    cm_s{n} = confusionmat(YPred{n}, YValidation{n});
    cm_sum = cm_s{n} + cm_sum;
end
%%
for c = 1:length(extracted_labels)
    precision(c) = cm_sum(c, c) / sum(cm_sum(:, c));
    recall(c) = cm_sum(c, c) / sum(cm_sum(c, :));
    f1Score(c) = 2 * (precision(c) * recall(c)) / (precision(c) + recall(c));

    fprintf('Class %s: Precision = %.4f, Recall = %.4f, F1 Score = %.4f\n', extracted_labels{c}, precision(c), recall(c), f1Score(c));
end
%%
labelsViruses = char(cellstr(extracted_labels)); 

cmChart1 = confusionchart(cm_sum,string(labelsViruses));
cmChart1.Normalization = 'column-normalized';
set(gca,'fontsize',16)

%%


classes ={'0' '1' '2' '3' '4' '6'};
objroc = rocmetrics(YValidation, scores, classes);
plot(objroc, ShowModelOperatingPoint=false, AverageROCType="micro")
objroc.AUC
objroc.Metrics
hl = legend;
hl.String{1} = 'Alpha';
hl.String{2} = 'Beta';
hl.String{3} = 'Delta';
hl.String{4} = 'Gamma';
hl.String{5} = 'Epsilon';
hl.String{6} = 'Iota';

