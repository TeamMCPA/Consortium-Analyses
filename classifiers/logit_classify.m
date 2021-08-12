function [classification, comparisons] = logit_classify(train_data, train_labels, test_data, test_labels, opts)
%% logistic regression wrapper
% takes in training data, training labels, testing data, and testing labels
% (unused) as well as opts struct with classification parameters

%% parse parameters
pairwise = opts.pairwise;
opts = rmfield(opts,'pairwise');
input = parse_opts(opts);

%% if pairwise classification
if pairwise
    
    % set up some parameters for our compairison loop
    number_classes = length(unique(train_labels));
    class_names = unique(train_labels);
    class_count = sum(strcmp(test_labels, class_names{1}));
    list_of_comparisons = combnk([1:number_classes],2);
    number_of_comparisons = size(list_of_comparisons,1);
    results_of_comparisons = cell((class_count*2), 2, number_of_comparisons);
    
    for this_comp = 1:number_of_comparisons
        
        % figure out what classes we need to keep for this round of
        % comparisons
        test_classes = list_of_comparisons(this_comp,:);
        test_class_names = class_names(test_classes);
        
        % select just the rows of the training data and labels that we will
        % need for this comparison
        train_dat = [train_data(strcmp(train_labels,test_class_names(1)),:);...
                     train_data(strcmp(train_labels,test_class_names(2)),:)];
        train_labs = [train_labels(strcmp(train_labels,test_class_names(1)));...
                      train_labels(strcmp(train_labels,test_class_names(2)))];
        
        
        % select just the rows of the test data and labels that we will
        % need for this comparison
        test_dat = [test_data(strcmp(test_labels,test_class_names(1)),:);...
                    test_data(strcmp(test_labels,test_class_names(2)),:)];
        test_labs = [test_labels(strcmp(test_labels,test_class_names(1)));...
                     test_labels(strcmp(test_labels,test_class_names(2)))];
        
        % create the model
        t = templateLinear('Learner', 'logistic', input{:});
        logit_model = fitcecoc(train_dat, train_labs, 'Learners', t);
        
        % perform the classification
        classification = predict(logit_model, test_dat);
        
        results_of_comparisons(:,1, this_comp) = classification;
        results_of_comparisons(:,2, this_comp) = test_labs;  
    end
    comparisons = list_of_comparisons;
    classification = results_of_comparisons;
else
    t = templateLinear('Learner', 'logistic', input{:});
    logit_model = fitcecoc(train_data, train_labels, 'Learners', t);
    
    classification = predict(logit_model, test_data);
    
    comparisons = nan(length(test_labels),1);
    unique_labels = unique(test_labels);
    for i = 1:length(unique_labels)
        comparisons(strcmp(test_labels, unique_labels{i})) = i;
    end
end


end
