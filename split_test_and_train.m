function [train_data, train_labels, test_data, test_labels] = split_test_and_train(fold_idx, conditions, pattern_data, event_types, final_dimensions, dimension_labels) 

% Set logical flags for indexing the conditions that will be compared.
% Loop through the whole list of conditions and create flags for each.
cond_flags = cell(length(conditions),1); % These are, for the moment, empty
%% step 1: find the locations of the conditions we would like to classify
for cond_idx = 1:length(conditions)
    if ischar(conditions{cond_idx}) || isstring(conditions{cond_idx}) || iscellstr(conditions{cond_idx})
        [~, ~, cond_flags{cond_idx}] = intersect(conditions{cond_idx},event_types);
    else
        cond_flags{cond_idx} = conditions{cond_idx};
    end
end

%% step 2: get the train data
% define indices of train data and test data
% find where the folding dimension is
group_vec = [1:size(pattern_data,ndims(pattern_data))];
group_vec = group_vec(group_vec~=fold_idx);

% first find the dimensions to concatenate
split_dims = cellfun(@(s) strsplit(s,'X'), final_dimensions, 'UniformOutput', false);
has_dims_to_concat = cellfun('size', split_dims,2) > 1;

if has_dims_to_concat
    concat_dims = split_dims{has_dims_to_concat};

    to_concat = [];
    for i = 1:length(concat_dims)
        idx = find(strcmp(dimension_labels, concat_dims{i}));
        to_concat = [to_concat, idx];
    end
    % get train data
    if ndims(pattern_data) == 3
        train_data = concatenate_dimensions(pattern_data([cond_flags{:}],:,group_vec), to_concat);
    else
        train_data = concatenate_dimensions(pattern_data([cond_flags{:}],:,:,group_vec), to_concat);
    end
else
    fprintf('this happened')
    if ndims(pattern_data) == 3
        train_data = pattern_data([cond_flags{:}],:,group_vec);
                size(train_data)

    else
        train_data = pattern_data([cond_flags{:}],:,:,group_vec);
        size(train_data)
    end
end

%% step 2: get test data

% subset out test
if ndims(pattern_data) == 3
    temp_test = pattern_data([cond_flags{:}],:,fold_idx);
else
    temp_test = pattern_data([cond_flags{:}],:,:,fold_idx);
end

if has_dims_to_concat
    % need to account for the assignment of 1 fold dropping the last dimension
    [max_value, max_ind] = max(to_concat);
    if max_value > ndims(temp_test)
        to_concat(max_ind) = [];
    end

    % create the test set
    test_data = concatenate_dimensions(temp_test, to_concat);
else
    test_data = temp_test;
end
    
%% step 4: get the class labels for test and train data

% first need to know how many times each label is represented in the data
train_label_repetitions = size(train_data, 1)/length(conditions);
test_label_repetitions = size(test_data, 1)/length(conditions);

% then create vector of labels
train_labels = repmat(conditions, 1,train_label_repetitions)';
test_labels = repmat(conditions, 1,test_label_repetitions)';
  
    
end
