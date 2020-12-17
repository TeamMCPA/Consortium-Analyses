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

%% Parse out the input data
input_struct = parse_inputs(MCP_struct, varargin{:});  

%% validate classification options
input_struct.opts_struct = validate_classification_options_input(input_struct, input_struct.suppress_warnings);

%% Setting up the combinations of feature subsets

% Create all possible subsets. If setsize is equal to the total number of
% features, there will only be one 'subset' which is the full feature
% array. If setsize is less than the total number of features, there will
% be n-choose-k subsets to analyze.
%
% The size of the subsets can grow extremely quickly with the size of
% incl_features. Consequently, there is a default max of 1000000 sets,
% which can be customized. If the total number of sets is larger than the
% max number of allowed sets, the list of sets will be subsampled.

% Determine how many sets will be generated. Can use this later for warning
% messages or other branching. Sets variable turns into a huge memory hog.

unmapped_sets = find_feature_sets(input_struct);
sets = map_features_to_sets(input_struct, unmapped_sets);


%% norm check - do we want to scale individual participant data?

if input_struct.scale_data
    MCP_struct = scale_individuals(MCP_struct, input_struct);
end

%% Build MCPA struct for all subjects in the MCP

mcpa_struct = MCP_to_MCPA(MCP_struct,...
     input_struct.incl_subjects,...
     input_struct.incl_features,...
     input_struct.incl_channels,...
     input_struct.time_window,...
     input_struct.baseline_window,...
     input_struct.hemoglobin);

% Subset patterns by session
inds = pad_dimensions(mcpa_struct.dimensions, 'session', input_struct.incl_sessions);
mcpa_struct.patterns = mcpa_struct.patterns(inds{:}); % Replace patterns matrix with the subsetted sessions data


%% summarize MCPA struct

% Step 2: Apply the desired function (e.g., @nanmean) for summarizing time
% window data. You can write custom functions to deal with time- and
% feature-domain data however you want. Default behavior is to apply the
% function along the first dimension of the MCPA pattern (instance) and then the second dimension (time),
% but this can also be changed.

%% first decide how we want to concatenate or average over our dimensions

% intermediary step: see if the user specified the summarizing dimensions. If not,
% recommend what dimensions to average over

if ~isempty(input_struct.summarize_dimensions) || ~isfield(input_struct, 'summarize_dimensions')
    summarize_dimensions = input_struct.summarize_dimensions;   
else
    isWithinSubjects = true;
    warning('summarize_dimensions not specified. Consulting recommend_dimensions.')
    
    [summarize_dimensions, ~] = recommend_dimensions(input_struct, isWithinSubjects);
    
    fprintf('Summarizing dimensions with %s:\n',func2str(input_struct.summary_handle))
    fprintf('%s ',summarize_dimensions{:})
    fprintf('\n')
end

% then see if the user specified the final dimensions the data should take
% before going into classification
if ~isempty(input_struct.final_dimensions) || ~isfield(input_struct, 'final_dimensions')
    final_dimensions = input_struct.final_dimensions;
else
    isWithinSubjects = true;
    warning('final_dimensions not specified. Consulting recommend_dimensions.')
    
    [~, final_dimensions] = recommend_dimensions(input_struct, isWithinSubjects);
   
    fprintf('The format the data will be in when it enters the classifier wrapper is: %s', final_dimensions{:});
    fprintf('\n')
end

%% then do the summarizing

if input_struct.verbose
    disp('Summarizing MCPA patterns with dimensions:');
    disp(strjoin(mcpa_struct.dimensions,' x '));
    disp(strjoin(cellfun(@num2str, num2cell(size(mcpa_struct.patterns)),'UniformOutput',false),' x '));
end

mcpa_summ = summarize_MCPA_Struct(input_struct.summary_handle,...
    mcpa_struct,...
    summarize_dimensions);

if input_struct.verbose
    disp('MCPA patterns have been summarized to:')
    disp(strjoin(mcpa_summ.dimensions,' x '));
    disp(strjoin(cellfun(@num2str, num2cell(size(mcpa_summ.patterns)),'UniformOutput',false),' x '));
end

%% Prep some basic parameters

n_subj = length(input_struct.incl_subjects);
n_sets = size(sets,1);
n_feature = length(input_struct.incl_features);
s = size(mcpa_summ.patterns);
n_sessions = s(end-1);
try n_cond = length(unique(input_struct.conditions)); catch, n_cond = length(input_struct.conditions); end

%% Set up the results structure which includes a copy of MCPA_pattern

allsubj_results = create_results_struct(false,...
    mcpa_summ,...
    input_struct,...
    sets,...
    n_subj,...
    n_sets,...
    n_feature,...
    n_cond,...
    final_dimensions);
                                          
stack = dbstack;
current_folding_function = stack.name;
allsubj_results.test_type = current_folding_function;
allsubj_results.approach = input_struct.approach;
allsubj_results.randomized = input_struct.randomized_or_notrand;

%% for one subject at a time...

for s_idx = 1:n_subj
    
    if input_struct.verbose
        fprintf('Running %g feature subsets for Subject %g / %g \n',n_sets,s_idx,n_subj);
    end
    tic;
    
    inds = pad_dimensions(allsubj_results.dimensions, 'subject', s_idx);
    subject_patterns = mcpa_summ.patterns(inds{:});
    
    %% K-Fold (kf)
    
    if strcmp(input_struct.approach, 'kf')   
        subject_labels = repmat(mcpa_summ.event_types, (size(subject_patterns,1)/8),1);

        % remove empty rows
        remove = [];
        for r = 1:size(subject_patterns,1)
            if sum(isnan(subject_patterns(r,:))) > 0
               remove = [remove; r];
            end
        end

        subject_patterns(remove,:) = [];
        subject_labels(remove) = [];

        % find balanced classes
        num_data = size(subject_patterns,1);
        num_in_fold = floor(num_data*input_struct.test_percent);
        num_folds = floor(num_data/num_in_fold);
        fold_end_idx_array = [num_in_fold: num_in_fold: num_in_fold*num_folds];
        fold_end_idx_array(end) = num_data;
        fold_start_idx_array = [1:num_in_fold: num_in_fold*num_folds];

        
        % should we randomize the rows?
        if strcmp(input_struct.randomized_or_notrand, 'randomized')
            rng('default')
            rand_inds = randperm(size(subject_patterns,1));
            subject_patterns(rand_inds, :);
            subject_labels(rand_inds);
        end
        
        
        % should we balance the classes?
        if input_struct.balance_classes
            num_repetitions = num_data/length(input_struct.conditions);
            new_kfold_mat = validate_balanced_classes(fold_start_idx_array,...
                fold_end_idx_array,...
                input_struct,...
                num_folds,...
                subject_labels,...
                num_repetitions);
        end
        
        % define where test and train labels are derived from
        event_types = subject_labels;

    %% Leave-One-Out (loo)
    elseif strcmp(input_struct.approach, 'loo')   
        % randomize trials, if needed
        if strcmp(input_struct.randomized_or_notrand, 'randomized') 
            new_index = randperm(size(subject_patterns,ndims(subject_patterns)));
            subject_patterns = subject_patterns(:,:,new_index);
        end
        
        % define where test and train labels are derived from
        event_types = mcpa_summ.event_types;
        
        % define the fold
        num_folds = length(MCP_struct(s_idx).Experiment.Runs);
    end
    
    %% Begin cross-validation
    for folding_idx = 1:num_folds         
        %% Define fold_idx: indices in subject_patterns that will be test data
        temp_set_results_cond = nan(n_cond,n_sets,n_feature);
        
        if strcmp(input_struct.approach, 'kf')
            fold = new_kfold_mat(folding_idx,~isnan(new_kfold_mat(folding_idx,:)));
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
        
        [train_data, train_labels, test_data, test_labels] = split_test_and_train(fold,...
            input_struct.conditions,...
            subject_patterns,... 
            event_types,...
            final_dimensions,...
            mcpa_summ.dimensions, [], []);
                
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
        
        if strcmp(input_struct.approach, 'kf') && (~isfield(input_struct, 'randomsubset')||isempty(input_struct.randomsubset))
            allsubj_results.kf_randomsubset = 'no: using the whole dataset';
            allsubj_results.kf_sampling = 'no sampling';  
            
        elseif strcmp(input_struct.approach, 'kf') && isfield(input_struct, 'randomsubset')
            n_test = size(test_data,1);
            n_train = size(train_data,1);
            
            subset_test = floor(input_struct.randomsubset*input_struct.test_percent);
            subset_train = floor(input_struct.randomsubset*(1-input_struct.test_percent));
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
        
        for set_idx = 1:min(n_sets,input_struct.max_sets)    

            if input_struct.verbose
                status_jump = floor(n_sets/20);
                if ~mod(set_idx,status_jump)
                    fprintf(' .')
                end
            end
            
            % Select the features for this subset
            set_features = sets(set_idx,:);

            %% Classify
            inds = pad_dimensions(final_dimensions, 'feature', set_features);
            [predicted_labels, comparisons] = input_struct.test_handle(...
                    train_data(inds{:}), ...
                    train_labels,...
                    test_data(inds{:}),...
                    test_labels,...
                    input_struct.opts_struct);

            %% Record results 

            if size(predicted_labels,2) > 1 % test labels will be a column vector if we don't do pairwise
                if s_idx==1 && set_idx == 1 && folding_idx == 1
                    if strcmp(input_struct.approach, 'kf')   
                        allsubj_results.accuracy_matrix = nan(n_cond,n_cond,min(n_sets,input_struct.max_sets),num_folds,n_subj);                             
                    elseif strcmp(input_struct.approach, 'loo')   
                        allsubj_results.accuracy_matrix = nan(n_cond,n_cond,min(n_sets,input_struct.max_sets),n_sessions,n_subj); 
                    end
                end

                if iscell(comparisons)
                    subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
                    comparisons = cellfun(@(x) find(strcmp(x,input_struct.conditions)),comparisons); 
                else
                    subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
                end

                % folding_idx should be the same size as however many folds you do 
                for comp = 1:size(comparisons,1)
                    if size(comparisons,2)==1
                        allsubj_results.accuracy_matrix(comparisons(comp,1),:,set_idx,folding_idx,s_idx) = subj_acc(comp);
                    else
                        allsubj_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),set_idx,folding_idx,s_idx) = subj_acc(comp);
                    end
                end
            else
                for cond_idx = 1:n_cond
                    temp_acc = cellfun(@strcmp,...
                    test_labels(strcmp(strjoin(string(input_struct.conditions{cond_idx}),'+'),test_labels)),... % known labels
                    predicted_labels(strcmp(strjoin(string(input_struct.conditions{cond_idx}),'+'),test_labels))...% classifier labels
                    );

                    temp_set_results_cond(cond_idx,set_idx,set_features) = nanmean(temp_acc);
                end
                for cond_idx = 1:n_cond
                    allsubj_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
                    allsubj_results.accuracy(cond_idx).subjXfeature(s_idx,:) = nanmean(temp_set_results_cond(cond_idx,:,:),2);
                    allsubj_results.accuracy(cond_idx).subjXsession(s_idx,folding_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
                end
            end

        end % end set_idx
        
    if input_struct.verbose
        fprintf(' %0.1f mins\n',toc/60);
    end
    end % end session loop
    
end % end subject loop


%end % end function
