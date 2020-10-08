function results_struct = classify_ParticipantLevel(results_struct)
%% nfold_classify_ParticipantLevel takes in a struct containing the MCPA summarized patterns 
% and classification parameters parsed from cross_validate,
% nest_cross_validate, or significance_test. It will then perform between
% subjects cross validation, where on each fold, one participant's data
% will be left out and the rest will used as training data for the
% classifier. It then returns the original struct with the accuracy of each
% fold. 

% created by Anna Herbolzheimer and Ben Zinszer summer 2020

%% Prep some basic parameters
n_subj = length(results_struct.incl_subjects);
n_sets = size(results_struct.subsets,1);
n_feature = length(results_struct.incl_features);
try n_cond = length(unique(results_struct.conditions)); catch, n_cond = length(results_struct.conditions); end

%% Folding & Dispatcher
for s_idx = 1:n_subj    
    if results_struct.verbose
        fprintf('Running %g feature subsets for Subject %g / %g',n_sets,s_idx,n_subj);
    end
    tic;
    
    %% Run over feature subsets
    temp_set_results_cond = nan(n_cond,n_sets,n_feature);
    
    %% Split data into test and train 
    % on each fold, one participant's data will be left out as the test
    % set, the rest of the participants data will be combined according to
    % final_dimensions. 

    [train_data, train_labels, test_data, test_labels] = split_test_and_train(s_idx,...
        results_struct.conditions,...
        results_struct.patterns,...
        results_struct.conditions,...
        results_struct.final_dimensions,...
        results_struct.dimensions, [], []);

    
    %% Run feature subsetting search
    for set_idx = 1:min(n_sets,results_struct.max_sets)
        %% Progress reporting bit (not important to function. just sanity)
        % Report at every 5% progress
        if results_struct.verbose
            status_jump = floor(n_sets/20);
            if ~mod(set_idx,status_jump)
                fprintf(' .')
            end
        end
        % Select the features for this subset
        set_features = results_struct.subsets(set_idx,:);
        
        %% Classify
        inds = subset_dimension(results_struct.final_dimensions, 'feature', set_features);
        [predicted_labels, comparisons] = results_struct.test_handle(...
                train_data(inds{:}), ...
                train_labels,...
                test_data(inds{:}),...
                test_labels,...
                results_struct.opts_struct);
            
        %% Record results 
        results_struct = save_accuracy(results_struct,...
            comparisons,...
            predicted_labels,...
            set_idx,...
            1,... 
            s_idx);

        
        %% Progress reporting
        if results_struct.verbose
            fprintf(' %0.1f mins\n',toc/60);
        end
    end % end set_idx loop
end % end subject loop
end