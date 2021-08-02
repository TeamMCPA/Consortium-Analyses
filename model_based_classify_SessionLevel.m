function allsubj_results = model_based_classify_SessionLevel(MCP_struct, semantic_model,semantic_model_labels, varargin)
%% function to classify sessions from individual participants
% n-fold within subjects cross-validation for n subjects with m sessions
% to classify individual sessions' average response patterns.


% This wrapper allows the user to choose if
% features will be averaged within-participants to produce a single
% participant-level observation or if individual events will be preserved.
% If the featuers are averaged within-participants, the training set
% is constrained to the number of participants minus 1. Otherwise the
% training set will be (participants - 1) * (number of instances of an
% object).

% Several parameters can be changed,
% including which functions are used to generate features and what
% classifier is trained. See Arguments below:
%
% Arguments:
% MCP_struct: either an MCP-formatted struct or the path to a Matlab file
% (.mat or .mcp) containing the MCP_struct.
% incl_features: features to include in the analysis. Default: all features
% incl_subjects: index of participants to include. Default: all participants
% baseline_window: [onset, offset] in seconds. Default [-5,0]
% time_window: [onset, offset] in seconds. Default [2,6]
% conditions: cell array of condition names / trigger #s. Default: {1,2}
% summary_handle: function handle (or char of function name) to specify how
% time-x-feature data should be summarized into features. Default: nanmean
% setsize: number of features to analyze (for subset analyses) Default: all
% test_handle: function handle for classifier. Default: mcpa_classify
% opts_struct: contains additional classifier options. Default: empty struct
% verbose: logical flag to report status updates and results. Default: true
% norm_data: option to perform some for of feature scaling. Defualt: false
% norm_withinSessions: if we do norm data, this allows the user to choose
% if we norm within individual sessions or across a participants session.
% It is only valid to norm across sessions for
% nfold_classify_ParticipantLevel. Default: true
% norm_function: function handle for how to scale the data. Default:
% minMax_scale
% minMax: what min and max to set the data to. Default: [0,1]
% summarize_dimensions: what dimensions and what order to summarize the
% dimensions of mcpa patterns by. Default behavior is to average dimensions
% in the order of this cell array. Default: {'instance', 'time'}


%% Load MCP struct if necessary
if isstring(MCP_struct) || ischar(MCP_struct)
    MCP_struct = load(MCP_struct,'-mat');
    varname = fieldnames(MCP_struct);
    MCP_struct = eval(['MCP_struct.' varname{1}]);
end

%% if data has already been summarized, leave as is. Otherwise, setup MCPA data and summarize it
if ~any(cellfun(@(x) strcmp(x, 'results_struct'), varargin(find(rem(1:length(varargin), 2)))))
    allsubj_results = setup_MCPA_data(MCP_struct,varargin);    
else
    allsubj_results = varargin{find(rem(1:length(varargin), 2))+1};
end   

%% Make sure that semantic model is a correlation matrix and if not, turn it into one
[semantic_model,semantic_model_labels] = validate_model_matrix(semantic_model, semantic_model_labels, allsubj_results.conditions, allsubj_results);

%% Prep some basic parameters
n_subj = length(allsubj_results.incl_subjects);
n_sets = size(allsubj_results.subsets,1);
n_feature = length(allsubj_results.incl_features);
s = size(allsubj_results.patterns);
n_sessions = s(end-1);
try n_cond = length(unique(allsubj_results.conditions)); catch, n_cond = length(allsubj_results.conditions); end

allsubj_results.approach = allsubj_results.approach;
allsubj_results.randomized = allsubj_results.randomized_or_notrand;

%% now begin the fold

for s_idx = 1:length(MCP_struct)
    
    if allsubj_results.verbose
        fprintf('Running %g feature subsets for Subject %g / %g \n',n_sets,s_idx,n_subj);
    end
    tic;
    
    % might want to get rid of this part later - purpose is to subset out
    % the participant for this round of CV
    if length(size(allsubj_results.patterns)) == 5
        subject_patterns = allsubj_results.patterns(:,:,:,:,s_idx);
    else
        subject_patterns = allsubj_results.patterns(:,:,:,s_idx);
    end
    
    %% do we have more than one session to classify on?
    % with current summarize dimensions, we can only proceed with
    % classification for 1 session if we are using KNN, SVM, or logistic
    % regression
    if n_sessions == 1
        %find something else for test train to be
        fold_dim = ndims(subject_patterns);
        num_data = size(subject_patterns,fold_dim);
        
        if ~isfield(allsubj_results, 'test_percent') || isempty(allsubj_results.test_percent)
            test_percent = .2;
        else
            test_percent = allsubj_results.test_percent;
        end
        num_in_fold = num_data*test_percent;
        num_folds = num_data/num_in_fold;
        one_session = true;
    else
        num_folds = length(MCP_struct(s_idx).Experiment.Runs);
        one_session = false;
    end
    
    for folding_idx = 1:num_folds
        %% Run over feature subsets
        temp_set_results_cond = nan(n_cond,n_sets,n_feature);
        
        if one_session
            fold_idx_end = folding_idx * num_in_fold;
            fold_idx_start = fold_idx_end - num_in_fold + 1;
            fold = fold_idx_start:fold_idx_end;
        else
            fold = folding_idx;
        end
                
         %% Folding & Dispatcher: Here's the important part
        % Right now, the data have to be treated differently for 2
        % conditions vs. many conditions. In MCPA this is because 2
        % conditions can only be compared in feature space (or, hopefully,
        % MNI space some day). If there are a sufficient number of
        % conditions (6ish or more), we abstract away from feature space
        % using RSA methods. Then classifier is trained/tested on the RSA
        % structures. This works for our previous MCPA studies, but might
        % not be appropriate for other classifiers (like SVM).
        
        [~, ~, test_data, test_labels] = split_test_and_train(fold,...
            allsubj_results.conditions,...
            subject_patterns,...
            allsubj_results.event_types,...
            allsubj_results.final_dimensions,...
            allsubj_results.dimensions, [], []);
        
        % permute the group labels if significance testing 
        if allsubj_results.permutation_test
            num_labels = length(semantic_model_labels);
            permuted_idx = randperm(num_labels)';
            semantic_model_labels = semantic_model_labels(permuted_idx);
        end 
                                                                            
        %% Run classifier and compare output with correct labels
        for set_idx = 1:min(n_sets,allsubj_results.max_sets)    
            %% Progress reporting bit (not important to function. just sanity)
            % Report at every 5% progress
            if allsubj_results.verbose
                status_jump = floor(n_sets/20);
                if ~mod(set_idx,status_jump)
                    fprintf(' .')
                end
            end
            % Select the features for this subset
            set_features = allsubj_results.subsets(set_idx,:);

             %% Classify
            inds = pad_dimensions(allsubj_results.final_dimensions, 'feature', set_features);
            [predicted_labels, comparisons] = allsubj_results.test_handle(...
                    semantic_model, ...
                    semantic_model_labels,...
                    test_data(inds{:}),...
                    test_labels,...
                    allsubj_results.opts_struct);

            %% Record results 
            if size(predicted_labels,2) > 1 % test labels will be a column vector if we don't do pairwise          
                subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
                nan_idx = cellfun(@(x) any(isnan(x)), predicted_labels(:,1,:), 'UniformOutput', false);
                subj_acc(:,:,[nan_idx{1,:,:}]) = nan;

                % Then loop through comparisons and save accuracy to the results struct
                for comp = 1:size(comparisons,1)
                    allsubj_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),set_idx, folding_idx, s_idx) = subj_acc(comp);
                end

            else
                subj_acc = double(strcmp(predicted_labels, test_labels));
                nan_idx = cellfun(@isnan, predicted_labels);
                subj_acc(nan_idx) = NaN;
                for cond_idx = 1:n_cond
                    cond_acc = nanmean(subj_acc(comparisons == cond_idx));
                    allsubj_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = cond_acc;
                    allsubj_results.accuracy(cond_idx).subjXfeature(s_idx,:) = cond_acc;
                    allsubj_results.accuracy(cond_idx).subjXsession(s_idx,folding_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
                end
            end


        end % end set_idx
    %% Progress reporting
    if allsubj_results.verbose
        fprintf(' %0.1f mins\n',toc/60);
    end
    end % end session loop
    
end % end subject loop


end % end function