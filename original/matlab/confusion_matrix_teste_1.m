%%
extracted_labels = cellfun(@(x) regexprep(x, '(\w+\s+\w+).*', '$1'), categories_samples, 'UniformOutput', false);

[predict_test, scores_teste] = classify(net,seqsCGR_teste);

% Get unique variant names and their indices
[categories_samples, ~, variant_indices] = unique(test_data.Variant);

% Assign labels based on the mapping
Labelteste_CGR = variant_indices;

% Convert to categorical
Labelteste_CGR = categorical(Labelteste_CGR);

cm_teste = confusionmat(predict_test, Labelteste_CGR);

labelsViruses_test = char(cellstr(extracted_labels)); 

test_chart = confusionchart(cm_teste, string(labelsViruses_test));

test_chart.Normalization = 'column-normalized';
set(gca, 'fontsize', 16)

%%
for c = 1:length(extracted_labels)
    precision(c) = cm_teste(c, c) / sum(cm_teste(:, c));
    recall(c) = cm_teste(c, c) / sum(cm_teste(c, :));
    f1Score(c) = 2 * (precision(c) * recall(c)) / (precision(c) + recall(c));

    fprintf('Class %s: Precision = %.4f, Recall = %.4f, F1 Score = %.4f\n', extracted_labels{c}, precision(c), recall(c), f1Score(c));
end

%%

classes ={'0' '1' '2' '3' '4' '6'};
objroc = rocmetrics(predict_test, scores_teste, classes);
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