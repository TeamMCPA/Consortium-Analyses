function results_struct = save_train_accuracy(results_struct, predicted_labels, correct_labels, set_idx, train_idx, s_idx)
%% Compare the labels output by the classifier to the known labels

if size(predicted_labels,2) > 1 % test labels will be a column vector if we don't do pairwise

    subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
    nan_idx = cellfun(@(x) any(isnan(x)), predicted_labels(:,1,:), 'UniformOutput', false);
    subj_acc(:,:,[nan_idx{1,:,:}]) = nan;

    % Then loop through comparisons and save accuracy to the results struct
    for comp = 1:size(correct_labels,1)
        results_struct.accuracy_matrix_train(correct_labels(comp,1),correct_labels(comp,2),set_idx,train_idx,s_idx) = subj_acc(comp);
    end

else
    subj_acc = strcmp(predicted_labels, correct_labels);

    for cond_idx = 1:n_cond
        cond_acc = nanmean(subj_acc(comparisons == cond_idx));
        results_struct.train_accuracy(cond_idx).subsetXsubj(:,s_idx) = cond_acc;
        results_struct.train_accuracy(cond_idx).subjXfeature(s_idx,:) = cond_acc;

    end
end


end