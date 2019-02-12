
function classification = mcpa_classify(model_data, model_labels, test_data, opts)

%% If the options struct is not provided, set default parameters
if ~exist('opts','var')
    opts = struct;
    opts.corr_stat = 'pearson';
end

%% Average across training data to get model features for each class
model_classes = unique(model_labels);
model_pattern1 = nanmean(model_data(strcmp(model_classes{1},model_labels),:),1)';
model_pattern2 = nanmean(model_data(strcmp(model_classes{2},model_labels),:),1)';

%% Perform correlation (default: Pearson) between two model patterns and the test patterns
% This is a quick and easy way to compute all the test items at once. 

% 1. Build a matrix with [ModelA, ModelB, Test1, Test2,... TestN]
% 2. Run correlation over all of them, and get a the matrix of corr coeffs.
% 3. First column is Model A vs. all others, Second column is Model B vs.
% all others.
corr_matrix = atanh(corr([model_pattern1,model_pattern2,test_data'],'type',opts.corr_stat,'rows','pairwise'));

test_model_corrs = corr_matrix(3:end,1:2);

%% Save out the classification results based on greatest correlation coefficient for each test pattern
% Initialize empty cell matrix for classifications
classification = cell(size(test_data,1),1);

% Classify based on greatest correlation coefficient
classification(test_model_corrs(:,1)>test_model_corrs(:,2)) = model_classes(1);
classification(test_model_corrs(:,1)<test_model_corrs(:,2)) = model_classes(2);

% If correlations to each model pattern are equal, enter NaN
classification(test_model_corrs(:,1)==test_model_corrs(:,2)) = {NaN};

end