function results_struct = save_accuracy(results_struct, correct_labels, predicted_labels, set_idx, session_idx, s_idx)
%% Compare the labels output by the classifier to the known labels

%% fill in accuracy matrix
if size(predicted_labels,2) > 1 % test labels will be a column vector if we don't do pairwise
    % Then see if the classification results are empty
    if ~iscell(predicted_labels(:,1,:)) && sum(isnan(predicted_labels(:,1,:)))==numel(predicted_labels(:,1,:))
        % if nans, then we'll need to save the subj_acc as nans
        if iscell(correct_labels) 
            % if correct_labels is cell array, we'll need to find the values in that array correspond to the condition location in event_types
            subj_acc = nan(length(predicted_labels(:,1,:)),1);
            correct_labels = cellfun(@(x) find(strcmp(x,mcpa_summ.event_types)),correct_labels);
        else
            subj_acc = nan(length(predicted_labels(:,1,:)),1);
        end
    elseif iscell(predicted_labels(:,1,:)) && sum(sum(cellfun(@(a) ~ischar(a), predicted_labels(:,1,:)))) == numel(predicted_labels(:,1,:))
        % if nans, then we'll need to save the subj_acc as nans
        if iscell(correct_labels) 
            % if correct_labels is cell array, we'll need to find the values in that array correspond to the condition location in event_types
            subj_acc = nan(length(predicted_labels(:,1,:)),1);
            correct_labels = cellfun(@(x) find(strcmp(x,mcpa_summ.event_types)),correct_labels);
        else
            subj_acc = nan(length(predicted_labels(:,1,:)),1);
        end
    else
        % if not nans, see how correct_labels is saved. 
        % If its a cell array, then we need to find where the values in that array correspond to the condition location in event_types
        % Then find the subject accuracy
        if iscell(correct_labels) 
            subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
            correct_labels = cellfun(@(x) find(strcmp(x,mcpa_summ.event_types)),correct_labels);
        else
            subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
        end
    end

    % Then loop through correct_labels and save accuracy to the results struct
    for comp = 1:size(correct_labels,1)
        if size(correct_labels,2)==1
            results_struct.accuracy_matrix(correct_labels(comp,1),:,set_idx,session_idx, s_idx) = subj_acc(comp);
        else
            results_struct.accuracy_matrix(correct_labels(comp,1),correct_labels(comp,2),set_idx,session_idx, s_idx) = subj_acc(comp);
        end
    end
else
    for cond_idx = 1:length(results_struct.conditions)
        temp_acc = cellfun(@strcmp,...
            correct_labels(strcmp(strjoin(string(p.Results.conditions{cond_idx}),'+'),correct_labels)),... % known labels
            predicted_labels(strcmp(strjoin(string(p.Results.conditions{cond_idx}),'+'),correct_labels))...% classifier labels
            );

        temp_set_results_cond(cond_idx,set_idx,set_features) = nanmean(temp_acc);
    end
    for cond_idx = 1:length(results_struct.conditions)
        results_struct.accuracy(cond_idx).subsetXsubj(:,s_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
        results_struct.accuracy(cond_idx).subjXfeature(s_idx,:) = nanmean(temp_set_results_cond(cond_idx,:,:),2);
        
        if isfield(allsubj_results.accuracy(cond_id), 'subjXsession')
            allsubj_results.accuracy(cond_idx).subjXsession(s_idx,session_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);           
        end
        
    end
end
end