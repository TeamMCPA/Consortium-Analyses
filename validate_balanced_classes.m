function new_kfold_mat = validate_balanced_classes(fold_start_idx_array, fold_end_idx_array, results_struct, num_folds, subject_labels)

%% start by getting the counts of each class in each fold
validation_mat = nan(length(results_struct.conditions), num_folds);
for i = 1:num_folds
    idx = fold_start_idx_array(i):fold_end_idx_array(i);
    validation_mat(:,i) = groupcounts(subject_labels(idx));
    kfold_mat(i,:) = idx; 
end

%% initialize a new matrix to add the fold indices to
new_kfold_mat = nan(size(kfold_mat,1), size(kfold_mat,2)+4);
new_kfold_mat(1:size(kfold_mat,1), 1:size(kfold_mat,2)) = kfold_mat;

%% how many instances of each class should there be?
class_num = mode(validation_mat, 'all');

%% move indices around to balance out the classes
rows_added_to = [];
for row = 1:length(results_struct.conditions)
    % see if the class needs to be balanced across folds
    current_set = validation_mat(row,:);                   
    x = find(current_set ~= class_num);
    
    if ~isempty(x)
        amount_to_even = abs(diff(current_set(x)))/length(x);
        which_to_skim = find(current_set==max(current_set));
        which_to_pad = find(current_set==min(current_set));

        class_to_move = results_struct.conditions{row};
        class_inds = find(strcmp(subject_labels(kfold_mat(which_to_skim,:)), class_to_move));
        ind_to_select = randperm(length(class_inds)); % randomly select an index to move

        dim =  new_kfold_mat(which_to_skim, (class_inds(ind_to_select(1:length(amount_to_even))))); % get the real value of the selected indices to move

        rows_added_to = [rows_added_to, which_to_pad]; % keep track of how much has been added to each row
        for v = 1:length(dim)
            % first pad
            where_to_add = sum(rows_added_to == which_to_pad(v));
            new_kfold_mat(which_to_pad,size(kfold_mat,2)+where_to_add) =  dim(v);

            % then skim
            x = find(new_kfold_mat(which_to_skim, :) == dim(v));
            new_kfold_mat(which_to_skim,x) = NaN;
        end
    end
end

end
