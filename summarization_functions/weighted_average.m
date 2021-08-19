function tmp_mat = weighted_average(similarity_matrices, training_mat, opts)
%% compute the weighted average for the (dis)similarity matrices
% input - similarity_matrices: (dis)similarity matrices for either the test
%                              or traing data (condition x condition x
%                              session x subject)
%       - training_mat: true/false flag for if this is the training data
%       - opts: classification options struct
% output averaged data in format: condition x condition x session x subject

    % get max number of trials
    total_trials = unique(sum(opts.trials_per_session,2,'omitnan'));
    if length(total_trials) > 1
        error('uneven number of trials. no solution for that yet')
    end
    
    % subset out the training and test data
    tmp_mat = nan(size(similarity_matrices));
    if training_mat
        model_trials_per_session = opts.trials_per_session(setdiff(1:length(opts.trials_per_session), opts.fold),:);
    else
        model_trials_per_session = opts.trials_per_session(opts.fold,:);
    end
    
    % compute the weighted average
    for subj = 1:size(similarity_matrices,4)
        for sess = 1:size(similarity_matrices,3)
            tmp_mat(:,:,sess,subj) = (model_trials_per_session(subj, sess)*similarity_matrices(:,:,sess, subj))/total_trials;           
        end
    end
    tmp_mat = sum(tmp_mat,3, 'omitnan');

end
