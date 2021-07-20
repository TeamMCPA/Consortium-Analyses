function [classification, comparisons] = rsa_classify(model_data, model_labels, test_data, test_labels, opts)
%% rsa_classify implements a correlation-based, similarity-space classifier
% following Zinszer, Bayet, Emberson, Raizada, & Aslin's (2018, Neurophotonics)
% method for all-possible-pairwise-comparisons and Zinszer, Anderson, Kang,
% Wheatley, & Raizada's (2016, JoCN) approach for n-way comparison
%
% opts: a struct that contains options for the classifier.
%
% opts.comparison_type determines how the features will be abstracted:
%    'comparison_type' if true, will use correlation-based similarity
%    space (larger values more similar). opts.metric can be set to 
%    'spearman', 'pearson', or 'kendall'.
%    'comparison_type' if false, will use pdist-based dissimilarity matrix
%    (larger values are less similar). Any of the distance metrics from the
%    pdist function (e.g., 'cosine', 'euclidean', 'seuclidean') can be set
%    for the opts.metric value.
%    Default setting is opts.similiarty_space=true, opts.metric='spearman'
% 
%    If similarity space is used, values are adjusted with Fisher's r-to-z
%    transform (the hyperbolic arctangent of the correlation coefficient).
%
% Only opts.exclusive=true mode is currently available.
%
% opts.pairwise if true will make all possible pairwise contrasts between
% the classes. If false, it will test all n-way comparison of class labels
% for permutation of labels with highest correlation to the training set.
%
% If opts.tiebreak is set to true (default), the correlation coefficients
% are randomly adjusted by <1% of the smallest observed difference to
% prevent exact matches and thus prevent ties in the classification. If
% opts.tiebreak is set to false, the classifier will prefer classes
% appearing earlier in the training set.
%
% If opts.verbose is set to true (default is false), Fisher-adjusted
% (hyperbolic arctangent) correlation matrices are displayed for all
% subjects, sessions, etc.


%% Pull a list of all the unique classes / conditions, preserving order
model_classes = unique(model_labels(:),'stable');

%% check to see the orders of test and train data - this is to see if they match
% get number of categories

[~,train_order] = sort(model_labels);
[~,test_order] = sort(test_labels);

% sort test data based on the re-ordered data
test_labs = test_labels(test_order);
test_dat = test_data(test_order, :,:,:);

model_labs = model_labels(train_order);
model_dat = model_data(train_order, :,:,:);
    
%% Build similarity structures
% Transform the model_data into similiarty or dissimilarity structures for each session by
% correlating between conditions or finding pairwise differences between conditions and then 
% average the similarity/dissimilarity structures together into a single model.

% if using a premade training model (like a semantic model), don't
% need to do anything to model_dat
% if we're doing leave one out (participant or session), first want
% to correlate the model data within session, then average across sessions

% build the matrices
if ndims(model_dat) == 2 && size(model_dat,1) == length(model_classes)
    model_mat = model_dat;
else
    model_mat = opts.similarity_function(model_dat, opts);
end

test_mat = opts.similarity_function(test_dat, opts);

% then average across sessions and subjects
if opts.weighted_average
    training_matrix = weighted_average(model_mat, true, opts);
    training_matrix = atanh(nanmean(training_matrix,4));
    
    % repeat for test matrix
    test_matrix = weighted_average(test_mat, false, opts);
    test_matrix = atanh(nanmean(test_matrix,4));
    
else
    % then average across each participants' sessions
    training_matrix = nanmean(model_mat,3);
    % then average across participants
    training_matrix = nanmean(training_matrix,4);
    training_matrix = atanh(training_matrix);
    
    % repeat for test matrix
    test_matrix = nanmean(test_mat,3);
    test_matrix = nanmean(test_matrix,4);
    test_matrix = atanh(test_matrix);
    
end

%% Visualize the matrices

if opts.verbose > 1
    rsa_visualize_similarity_matrices(model_data, test_data, model_mat, test_mat)    
end

if opts.verbose
    rsa_visualize_test_and_training_matrices(training_matrix, model_labels, test_matrix, test_labels)
end

%% Sanity Check
if sum(isnan(test_matrix(:)))==(numel(test_matrix)-8) || sum(isnan(test_matrix(:)))==numel(test_matrix) || sum(isnan(training_matrix(:)))==numel(training_matrix)
     warning('One or both input matrices contains all NaN values. this one triggered');
end

%% Save out the classification results based on greatest correlation coefficient for each test pattern
% Initialize empty cell matrix for classifications

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
        results_of_comparisons(perm_idx) = corr(train_vec,tmp_test_vec,'rows','pairwise');
    end
    
    % Choose the best of all the permutations of labels
    [rating, best_perm] = max(results_of_comparisons);
    classification = model_classes(list_of_comparisons(best_perm,:));
    
    %accuracy = strcmp(classification,test_labels);
    comparisons = test_labels';
    
    % put the labels back in the order they were put in as
    [~,reorder_test] = sort(test_order);
    classification = classification(reorder_test);
    comparisons = reorder_test;
    
    
else
    
    [accuracy, comparisons] = pairwise_rsa_test(test_matrix,training_matrix);
    
    classification = comparisons;
    
    for acc = 1:length(accuracy)
        if isnan(accuracy(acc))
            classification(acc,:) = nan;
        else
            if ~accuracy(acc)
                classification(acc,:) = classification(acc,end:-1:1);
            end
        end
    end
   
    classification = num2cell(classification);
    for class = 1:size(classification,1)
        if ~any(isnan([classification{class,:}]))
            classification(class,:) = model_classes([classification{class,:}]);
        end
    end
    
    true_labels = model_classes(comparisons); 
    results_of_comparisons = cell(size(classification,2), 2, size(classification,1));
    for comp = 1:size(classification,1)
        results_of_comparisons(:,1,comp) = classification(comp,:)';
        results_of_comparisons(:,2,comp) = true_labels(comp,:)';
    end
    classification = results_of_comparisons;
    
    
end


end