function [classification rating] = rsa_classify(model_data, model_labels, test_data, test_labels, opts)
%% rsa_classify implements a correlation-based, similarity-space classifier
% following Zinszer, Bayet, Emberson, Raizada, & Aslin's (2018, Neurophotonics)
% method for all-possible-pairwise-comparisons and Zinszer, Anderson, Kang,
% Wheatley, & Raizada's (2016) approach for n-way comparison
%
% Correlation statistic (pearson, spearman, or kendall) may be selected
% using opts.corr_stat, e.g., default is opts.corr_stat='spearman'
%
% Only 'exclusive' mode is available.
%
% If opts.tiebreak is set to true (default), the correlation coefficients
% are adjusted by <1% of the smallest observed difference to prevent exact
% matches and thus prevent ties in the classification. If opts.tiebreak is
% set to false, the classifier will prefer classes appearing earlier in the
% training set.
%
% All-possible-pairwise comparison (opts.pairwise) is under development.

%% If the options struct is not provided, set default parameters
if ~exist('opts','var') || isempty(opts)
    opts = struct;
    opts.corr_stat = 'spearman';
    opts.exclusive = true;
    opts.pairwise = false;
    opts.tiebreak = true;
end

model_classes = unique(model_labels(:),'stable');

%% Build similarity structures
% Transform the model_data into similiarty structures for each session by
% correlating between conditions and then average the similarity structures
% together into a single model.

if length(size(model_data))>4
    % here assuming model_data are time-x-cond-x-chan-x-instance-x-session-x-subject
    model_data = squeeze(mean(model_data,4));
    
    % now time-x-cond-x-chan-x-session-x-subject
    model_data = squeeze(mean(model_data,1));
    
    % now cond-x-chan-x-session-x-subject
    model_correl = nan(size(model_data,2),size(model_data,2),size(model_data,3),size(model_data,4));
    for sub_idx = 1:size(model_data,4)
        for ses_idx = 1:size(model_data,3)
            model_correl = atanh(corr(model_data(:,:,ses_idx,sub_idx),'rows','pairwise','type','Spearman'));
        end
    end
    
    % now cond-x-cond-x-session-x-subject
    training_matrix = nanmean(nanmean(model_correl,3),4);
end

% Transform the test_data into similiarty structures for each session by
% correlating between conditions and then average the similarity structures
% together into a single structure.
if length(size(test_data))>4
    % here assuming model_data are time-x-cond-x-chan-x-instance-x-session-x-subject
    test_data = squeeze(mean(test_data,4));
    
    % now time-x-cond-x-chan-x-session-x-subject
    test_data = squeeze(mean(test_data,1));
    
    % now cond-x-chan-x-session-x-subject
    test_correl = nan(size(test_data,2),size(test_data,2),size(test_data,3),size(test_data,4));
    for sub_idx = 1:size(test_data,4)
        for ses_idx = 1:size(test_data,3)
            test_correl = atanh(corr(test_data(:,:,ses_idx,sub_idx),'rows','pairwise','type','Spearman'));
        end
    end
    
    % now cond-x-cond-x-session-x-subject
    test_matrix = nanmean(nanmean(test_correl,3),4);
end

%% Sanity Check
if sum(isnan(test_matrix(:)))==numel(test_matrix) || sum(isnan(training_matrix(:)))==numel(training_matrix)
    error('One or both input matrices contains all NaN values. I quit!');
end

%% Save out the classification results based on greatest correlation coefficient for each test pattern
% Initialize empty cell matrix for classifications
classification = cell(size(test_data,1),1);

if ~isfield(opts,'pairwise') || ~opts.pairwise
    
    % Number of stimulus classes to compare, based on the number of rows
    number_classes = size(test_matrix,1);
    
    % Too many permutations if more than 10 classes. Throw error and leave.
    if number_classes>10, error('n-way classification does not work for more than 10 classes. Use opts.pairwise=true'); end
    
    % Generate a list of every pairwise comparison and the results array
    list_of_comparisons = perms(1:number_classes);
    number_of_comparisons = size(list_of_comparisons,1);
    results_of_comparisons = nan(number_of_comparisons,1);
    
    % Convert training set to a vector for the comparisons
    train_vec = training_matrix(logical(tril(ones(size(training_matrix)),-1)));
    
    % Iterate through the different Permutations of test
    for perm_idx = 1:number_of_comparisons
        tmp_test_matrix = test_matrix(list_of_comparisons(perm_idx,:),list_of_comparisons(perm_idx,:));
        tmp_test_vec = tmp_test_matrix(logical(tril(ones(size(tmp_test_matrix)),-1)));
        results_of_comparisons(perm_idx) = corr(atanh(train_vec),atanh(tmp_test_vec),'rows','pairwise');
    end

    % Choose the best of all the permutations of labels
    [rating, best_perm] = max(results_of_comparisons);
    classification = model_classes(list_of_comparisons(best_perm,:));
    
else
    [accuracy, comparisons] = pairwise_rsa_test(test_matrix,training_matrix);
    
    %to-do: figure out how to create output-able data here.
end
