function allsubj_results = nfold_classify_WithinSubjects(MCP_struct, varargin)
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
%   (.mat or .mcp) containing the MCP_struct.
% incl_features: features to include in the analysis. Default: all features
% incl_subjects: index of participants to include. Default: all participants
% baseline_window: [onset, offset] in seconds. Default [-5,0]
% time_window: [onset, offset] in seconds. Default [2,6]
% conditions: cell array of condition names / trigger #s. Default: {1,2}
% summary_handle: function handle (or char of function name) to specify how
%   time-x-feature data should be summarized into features. Default: nanmean
% setsize: number of features to analyze (for subset analyses) Default: all
% test_handle: function handle for classifier. Default: mcpa_classify
% opts_struct: contains additional classifier options. Default: empty struct
% verbose: logical flag to report status updates and results. Default: true
% scale_data: option to perform some for of feature scaling. Defualt: false
% scale_withinSessions: if we do norm data, this allows the user to choose
%   if we norm within individual sessions or across a participants session.
%   It is only valid to norm across sessions for
%   nfold_classify_ParticipantLevel. Default: true
% scale_function: function handle for how to scale the data. Default:
%   minMax_scale
% minMax: what min and max to set the data to. Default: [0,1]
% summarize_dimensions: what dimensions and what order to summarize the
%   dimensions of mcpa patterns by. Default behavior is to average dimensions
%   in the order of this cell array. Default: {'instance', 'time'}
% appraoch: "loo" for leave-one-out; "kf" for k-fold
% randomized_or_notrand: Select whether or not to randomize the sessions 
%   (for loo) or trials (sessions x repetitions; for kf). Default: 'notrand'
% test_percent: if kf, percentage of data used for testing. Default = 0.2
% randomsubset: if randomly subset ths # of observations. Default = [], 
%   which uses all data 


% Created by: Anna Herbolzheimer and Ben Zinszer 2019
% edited by: Sori Baek and Anna Herbolzheimer 2020


%% Load MCP struct if necessary

if isstring(MCP_struct) || ischar(MCP_struct)
    MCP_struct = load(MCP_struct,'-mat');
    varname = fieldnames(MCP_struct);
    MCP_struct = eval(['MCP_struct.' varname{1}]);
end

%% if data has already been summarized, leave as is. Otherwise, setup MCPA data and summarize it
if ~any(cellfun(@(x) strcmp(x, 'results_struct'), varargin(find(rem(1:length(varargin), 2)))))
    allsubj_results = varargin{2};
    
else
    allsubj_results = setup_MCPA_data(MCP_struct,varargin);
end

%% Prep some basic parameters

n_subj = length(allsubj_results.incl_subjects);
n_sets = size(allsubj_results.subsets,1);
n_feature = length(allsubj_results.incl_features);
s = size(allsubj_results.patterns);
n_sessions = s(end-1);
try n_cond = length(unique(allsubj_results.conditions)); catch, n_cond = length(allsubj_results.conditions); end

%% for one subject at a time...
for s_idx = 1:n_subj
    
    if allsubj_results.verbose
        fprintf('Running %g feature subsets for Subject %g / %g \n',n_sets,s_idx,n_subj);
    end
    tic;
    
    inds = pad_dimensions(allsubj_results.dimensions, 'subject', s_idx);
    subject_patterns = allsubj_results.patterns(inds{:});
    
    %% set up indexing for either k-fold(kf) or leave one out (loo) CV
    %% K-Fold (kf)
    if strcmp(allsubj_results.approach, 'kf')   
        subject_labels = repmat(allsubj_results.event_types, (size(subject_patterns,1)/8),1);

        % remove empty rows
        remove = [];
        for r = 1:size(subject_patterns,1)
            if sum(isnan(subject_patterns(r,:))) > 0
               remove = [remove; r];
            end
        end

        subject_patterns(remove,:) = [];
        subject_labels(remove) = [];

        % find indices for each fold
        num_data = size(subject_patterns,1);
        num_in_fold = floor(num_data*allsubj_results.test_percent);
        num_folds = floor(num_data/num_in_fold);
        fold_end_idx_array = [num_in_fold: num_in_fold: num_in_fold*num_folds];
        fold_end_idx_array(end) = num_data;
        fold_start_idx_array = [1:num_in_fold: num_in_fold*num_folds];

        
        % should we randomize the rows?
        if strcmp(allsubj_results.randomized_or_notrand, 'randomized')
            rng('default')
            rand_inds = randperm(size(subject_patterns,1));
            subject_patterns(rand_inds, :);
            subject_labels(rand_inds);
        end
        
        
        % should we balance the classes?
        if allsubj_results.balance_classes
            num_repetitions = num_data/length(allsubj_results.conditions);
            new_kfold_mat = validate_balanced_classes(fold_start_idx_array,...
                fold_end_idx_array,...
                allsubj_results,...
                num_folds,...
                subject_labels,...
                num_repetitions);
        end
        
        % define where test and train labels are derived from
        event_types = subject_labels;

    %% Leave-One-Out (loo)
    elseif strcmp(allsubj_results.approach, 'loo')   
        % randomize trials, if needed
        if strcmp(allsubj_results.randomized_or_notrand, 'randomized') 
            new_index = randperm(size(subject_patterns,ndims(subject_patterns)));
            subject_patterns = subject_patterns(:,:,new_index);
        end
        
        % define where test and train labels are derived from
        event_types = allsubj_results.event_types;
        
        % define the fold
        num_folds = size(subject_patterns, ndims(subject_patterns));
    end
    
    %% Begin cross-validation
    for folding_idx = 1:num_folds         
        %% Define fold_idx: indices in subject_patterns that will be test data
        temp_set_results_cond = nan(n_cond,n_sets,n_feature);
        
        if strcmp(allsubj_results.approach, 'kf')
            fold = new_kfold_mat(folding_idx,~isnan(new_kfold_mat(folding_idx,:)));
        else
            fold = folding_idx;
        end

        %% Split data into test and train 
        % on each fold, one participant's data will be left out as the test
        % set, the rest of the participants data will be combined according to
        % final_dimensions. 
        
        [train_data, train_labels, test_data, test_labels] = split_test_and_train(fold,...
            allsubj_results.conditions,...
            subject_patterns,... 
            event_types,...
            {},...
            allsubj_results.dimensions, [], []);
                
        % Skip an iteration if all of train_data or all of test_data is NaN
        if sum(~isnan(train_data(:)))==0 || sum(~isnan(test_data(:)))==0
            fprintf('skip sub %d fold %d \n', s_idx, folding_idx);
            continue; 
        end
                   
        %% Take a random subset, if needed
        
        % Only for K-fold. If what you're trying to sample is greater
        % than the available dataset, this function will sample with
        % replacement within training data and within testing data separately,
        % so that no single data point is used in both training and testing
        % datasets. If what you'er trying to sample is less than the
        % available dataset, this function will sample without replacement.
        
        if strcmp(allsubj_results.approach, 'kf') && (~isfield(allsubj_results, 'randomsubset')||isempty(allsubj_results.randomsubset))
            allsubj_results.kf_randomsubset = 'no: using the whole dataset';
            allsubj_results.kf_sampling = 'no sampling';  
            
        elseif strcmp(allsubj_results.approach, 'kf') && isfield(allsubj_results, 'randomsubset')
            n_test = size(test_data,1);
            n_train = size(train_data,1);
            
            subset_test = floor(allsubj_results.randomsubset*allsubj_results.test_percent);
            subset_train = floor(allsubj_results.randomsubset*(1-allsubj_results.test_percent));
            allsubj_results.kf_randomsubset = sprintf('yes: %d out of %d test data and %d out of %d train data', subset_test, n_test, subset_train, n_train); 

            if subset_test <= n_test
                sampled_test_data = randsample(1:n_test,subset_test,false);
                sampled_train_data = randsample(1:n_train,subset_train,false);
                
                test_data = test_data(sampled_test_data,:);
                test_labels = test_labels(sampled_test_data);
                train_data = train_data(sampled_train_data,:);
                train_labels = train_labels(sampled_train_data);
                allsubj_results.kf_sampling = 'sample without replacement';  
                
            elseif  subset_test > n_test 
                sampled_test_data = randsample(1:n_test,subset_test,true);
                sampled_train_data = randsample(1:n_train,subset_train,true);
                
                test_data = test_data(sampled_test_data,:);
                test_labels = test_labels(sampled_test_data);
                train_data = train_data(sampled_train_data,:);
                train_labels = train_labels(sampled_train_data);
                
                allsubj_results.kf_sampling = 'sample with replacement';   
            end
            
        end

        %% Run classifier and compare output with correct labels        
        for set_idx = 1:min(n_sets,allsubj_results.max_sets)    

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
                    train_data(inds{:}), ...
                    train_labels,...
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
        
    if allsubj_results.verbose
        fprintf(' %0.1f mins\n',toc/60);
    end
    end % end session loop
    
end % end subject loop


%end % end function
