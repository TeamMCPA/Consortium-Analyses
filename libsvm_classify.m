function [classification, comparisons] = libsvm_classify(train_data, train_labels, test_data, test_labels, opts)
%% libsvm_classify implements a linear SVM classifier using the libSVM library
%
% Parameters:
% opts.pairwise = true, false (default = multiclass)
% opt.libsvm_options = a string of options in the same format as that of
% libSVM (default: '-s 0 -t 0 -q' i.e. -s 0 C-SVM, -t 0 linear, -q no output, with libSVM default values otherwise)
% no "spellcheck" of those options is made, and the user should ensure that
% they choose those options in a way that is appropriate for their
% analytic goals.

%% If the options struct is not provided, set default parameters
%% determine if we use pairwise comparisons
% pairwise option seemed to be forced to be true somewhere before this fucntion is called
% using manual pairwise value for debugging multiclass
if ~exist('opts','var') || ~isfield(opts, 'pairwise') || isempty(opts)
    pairwise = false;
else
    pairwise = opts.pairwise;
    opts = rmfield(opts,'pairwise');
end

if ~exist('pairwise','var')
    pairwise = false;
end

%libsvm options
if ~exist('opts','var') ||~isfield(opts, 'libsvm_options')|| isempty(opts)
    opts.libsvm_options = '-s 0 -t 0 -q';
end

%% parse out the classification parameters
input = parse_opts(opts);

%% Path
a=dir("libsvm*");
addpath(genpath(a.name));

%% remove NaN data in training data (libSVM can't deal with it)
nan_rows = sum(isnan(train_data),2)>0;
if sum(nan_rows)>0
    train_data = train_data(~nan_rows,:);
    train_labels = train_labels(~nan_rows);
    warning([num2str(sum(nan_rows)), 'rows with empty data out of ', num2str(length(nan_rows)) ,' were deleted in the training data']);
end
nan_cols = sum(isnan(train_data),1)>0;
if sum(nan_cols)>0
    train_data = train_data(:,~nan_cols);
    warning([num2str(sum(nan_cols)), 'columns with empty data out of ', num2str(length(nan_cols)) ,' were deleted in the training data']);
end
% same in test data
nan_rows = sum(isnan(test_data),2)>0;
if sum(nan_rows)>0
    test_data = test_data(~nan_rows,:);
    test_labels = test_labels(~nan_rows);
    warning([num2str(sum(nan_rows)), 'rows with empty data out of ', num2str(length(nan_rows)) ,' were deleted in the test data']);
end
nan_cols = sum(isnan(test_data),1)>0;
if sum(nan_cols)>0
    test_data = test_data(:,~nan_cols);
    warning([num2str(sum(nan_cols)), 'columns with empty data out of ', num2str(length(nan_cols)) ,' were deleted in the test data']);
end

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
        
        if isempty(train_dat) || isempty(test_dat)
            warning('Empty training or test data! Skipping...')
        else
            
        %% Convert labels to integer
        train_labs_int = 1+strcmp(train_labs,test_class_names{2});%pairwise 1 vs 2
        test_labs_int  = 1+strcmp(test_labs,test_class_names{2});%pairwise 1 vs 2
        
        %% Classifier training goes here
        model = svmtrain(train_labs_int, train_dat, opts.libsvm_options);  %#ok<SVMTRAIN> 

        %% Classifier testing goes here
        classification_int = svmpredict(test_labs_int, test_dat, model);%pairwise 1 vs 2
        classification = test_class_names(classification_int); %pairwise using label names
        
        % Save the predicted labels (`classification`) to first column
        results_of_comparisons(1:length(classification),1, this_comp) = classification;
        % Save the true test labels (`test_labs`) to second column
        results_of_comparisons(1:length(test_labs),2, this_comp) = test_labs;
        end
    end
    comparisons = list_of_comparisons;
    classification = results_of_comparisons;
    
%% if not pairwise classification
else

    %% multiclass (using libSVM default method)
    class_names = unique(train_labels);
    train_labels_int=NaN(size(train_labels));
    test_labels_int=NaN(size(test_labels));
    
    for i = 1:length(class_names)
         train_labels_int(strcmp(train_labels,class_names{i}))=i;
         test_labels_int(strcmp(test_labels,class_names{i}))=i;
    end
        
    if isempty(train_data) || isempty(test_data)
        classification=[];comparisons=[];
        warning('Empty training or test data! Skipping...')
    else
        %% Classifier training goes here
        model = svmtrain(train_labels_int, train_data, opts.libsvm_options); %#ok<SVMTRAIN>
        %% Classifier testing goes here
        classification_int = svmpredict(test_labels_int, test_data, model);
        classification = class_names(classification_int); %back to str
        
        comparisons = nan(length(test_labels),1);
        unique_labels = unique(test_labels);
        for i = 1:length(unique_labels)
            comparisons(strcmp(test_labels, unique_labels{i})) = i;
        end
    end
end

end
