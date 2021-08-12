function allsubj_results = model_based_classify_ParticipantLevel(MCP_struct, semantic_model, semantic_model_labels, varargin)
%% model_based_classify_ParticipantLevel takes an MCP struct and performs
% RSA classification for n subjects to classify individual
% participants' average response patterns using a semantic model. This wrapper assumes that
% features will be averaged within-participants to produce a single
% participant-level observation. Several parameters can be changed,
% including which functions are used to generate features. See Arguments below:
%
% Arguments:
% MCP_struct: either an MCP-formatted struct or the path to a Matlab file
% (.mat or .mcp) containing the MCP_struct.
% semantic_model: The model we want to use for classification
% incl_channels: channels to include in the analysis. Default: all channels
% incl_subjects: index of participants to include. Default: all participants
% baseline_window: [onset, offset] in seconds. Default [-5,0]
% time_window: [onset, offset] in seconds. Default [2,6]
% conditions: cell array of condition names / trigger #s. Default: {1,2}
% summary_handle: function handle (or char of function name) to specify how
% time-x-channel data should be summarized into features. Default: nanmean
% setsize: number of channels to analyze (for subset analyses) Default: all
% test_handle: function handle for classifier. Default: mcpa_classify
% opts_struct: contains additional classifier options. Default: empty struct
% verbose: logical flag to report status updates and results. Default: true

% dimensions for the accuracy matrix in results struct is cond x cond x
% channel subset x subject

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
try n_cond = length(unique(allsubj_results.conditions)); catch, n_cond = length(allsubj_results.conditions); end


%% Begin the n-fold process: Select one test subj at a time from MCPA struct
for s_idx = 1:length(allsubj_results.incl_subjects)
    if allsubj_results.verbose
        fprintf('Running %g feature subsets for Subject %g / %g',n_sets,s_idx,n_subj);
    end
    tic;
    
    %% Run over feature subsets
    temp_set_results_cond = nan(n_cond,n_sets,n_feature);
    
    %% Folding & Dispatcher: Here's the important part
    % Right now, the data have to be treated differently for 2
    % conditions vs. many conditions. In MCPA this is because 2
    % conditions can only be compared in feature space (or, hopefully,
    % MNI space some day). If there are a sufficient number of
    % conditions (6ish or more), we abstract away from feature space
    % using RSA methods. Then classifier is trained/tested on the RSA
    % structures. This works for our previous MCPA studies, but might
    % not be appropriate for other classifiers (like SVM).
       
    [~, ~, test_data, test_labels] = split_test_and_train(s_idx,...
        allsubj_results.conditions,...
        allsubj_results.patterns,...
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
                allsubj_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),set_idx, s_idx) = subj_acc(comp);
            end

        else
            subj_acc = double(strcmp(predicted_labels, test_labels));
            nan_idx = cellfun(@isnan, predicted_labels);
            subj_acc(nan_idx) = NaN;
            for cond_idx = 1:n_cond
                cond_acc = nanmean(subj_acc(strcmp(comparisons, test_labels(cond_idx))));
                allsubj_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = cond_acc;
                allsubj_results.accuracy(cond_idx).subjXfeature(s_idx,:) = cond_acc;
            end
        end
                  
        %% Progress reporting
        if allsubj_results.verbose
            fprintf(' %0.1f mins\n',toc/60);
        end
    end % end set_idx loop
end % end subject loop

end
