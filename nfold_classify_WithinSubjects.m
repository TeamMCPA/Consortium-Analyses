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
% n_randomsubset: if randomly subset ths # of observations. Default = [], 
%   which uses all data 


%% Load MCP struct if necessary

if isstring(MCP_struct) || ischar(MCP_struct)
    MCP_struct = load(MCP_struct,'-mat');
    varname = fieldnames(MCP_struct);
    MCP_struct = eval(['MCP_struct.' varname{1}]);
end

%% Parse out the input data

p = parse_inputs(MCP_struct, varargin{:}); % parse_inputs would have a new variable indicating what value of k we have for later

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

unmapped_sets = find_feature_sets(p.Results);
sets = map_features_to_sets(p, unmapped_sets);


%% norm check - do we want to scale individual participant data?

if p.Results.scale_data
    MCP_struct = scale_individuals(MCP_struct, p.Results);
end

%% Build MCPA struct for all subjects in the MCP

mcpa_struct = MCP_to_MCPA(MCP_struct,...
     p.Results.incl_subjects,...
     p.Results.incl_features,...
     p.Results.incl_channels,...
     p.Results.time_window,...
     p.Results.baseline_window,...
     p.Results.oxy_or_deoxy);

% Subset patterns by session
inds = repmat({':'},1,ndims(mcpa_struct.patterns)); % Create index structure with all-elements in all-dimensions
inds{strcmp(mcpa_struct.dimensions,'session')} = p.Results.incl_sessions; % In whichever dimension matches 'session', substitute the incl_sessions vector
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

if ~isempty(p.Results.summarize_dimensions) || ~isfield(p.Results, 'summarize_dimensions')
    summarize_dimensions = p.Results.summarize_dimensions;
    
elseif strcmp(p.Results.approach, 'kf')
    summarize_dimensions = {'time'};
    
else
    isWithinSubjects = true;
    warning('summarize_dimensions not specified. Consulting recommend_dimensions.')
    
    [summarize_dimensions, ~] = recommend_dimensions(p.Results, isWithinSubjects);
    
    fprintf('Summarizing dimensions with %s:\n',func2str(p.Results.summary_handle))
    fprintf('%s ',summarize_dimensions{:})
    fprintf('\n')
end

% then see if the user specified the final dimensions the data should take
% before going into classification
if ~isempty(p.Results.final_dimensions) || ~isfield(p.Results, 'final_dimensions')
    final_dimensions = p.Results.final_dimensions;
else
    isWithinSubjects = true;
    warning('final_dimensions not specified. Consulting recommend_dimensions.')
    
    [~, final_dimensions] = recommend_dimensions(p.Results, isWithinSubjects);
   
    fprintf('The format the data will be in when it enters the classifier wrapper is: %s', final_dimensions{:});
    fprintf('\n')
end

%% then do the summarizing

if p.Results.verbose
    disp('Summarizing MCPA patterns with dimensions:');
    disp(strjoin(mcpa_struct.dimensions,' x '));
    disp(strjoin(cellfun(@num2str, num2cell(size(mcpa_struct.patterns)),'UniformOutput',false),' x '));
end

mcpa_summ = summarize_MCPA_Struct(p.Results.summary_handle,...
    mcpa_struct,...
    summarize_dimensions);

if p.Results.verbose
    disp('MCPA patterns have been summarized to:')
    disp(strjoin(mcpa_summ.dimensions,' x '));
    disp(strjoin(cellfun(@num2str, num2cell(size(mcpa_summ.patterns)),'UniformOutput',false),' x '));
end

%% Prep some basic parameters

n_subj = length(p.Results.incl_subjects);
n_sets = size(sets,1);
n_feature = length(p.Results.incl_features);
s = size(mcpa_summ.patterns);
n_sessions = s(end-1);
try n_cond = length(unique(p.Results.conditions)); catch, n_cond = length(p.Results.conditions); end

%% Set up the results structure which includes a copy of MCPA_pattern

allsubj_results = create_results_struct(false,...
    mcpa_summ,...
    p,...
    sets,...
    n_subj,...
    n_sets,...
    n_feature,...
    n_cond,...
    final_dimensions);
                                          
stack = dbstack;
current_folding_function = stack.name;
allsubj_results.test_type = current_folding_function;

%% for one subject at a time...

for s_idx = 1:n_subj
    
    if p.Results.verbose
        fprintf('Running %g feature subsets for Subject %g / %g \n',n_sets,s_idx,n_subj);
    end
    tic;
    
    %% K-Fold (kf)
    
    if strcmp(p.Results.approach, 'kf')   
        allsubj_results.approach = 'k-fold (kf)';
        if length(size(mcpa_summ.patterns)) == 5
            subject_patterns = mcpa_summ.patterns(:,:,:,:,s_idx);
        else
            subject_patterns = mcpa_summ.patterns(:,:,:,s_idx);
        end % take only data from one subject 
        
        % concatenate all sessions & repetitions for each subject's condition and channel
        if n_sessions > 1 
            subject_patterns = concatenate_dimensions(subject_patterns, [ndims(subject_patterns)-1,ndims(subject_patterns)]); % concatenate all sessions together
        end
        
        allsubj_results.dimensions = {'condition', 'feature', 'repetitionXsession'};
        mcpa_summ.dimensions = {'condition', 'feature', 'repetitionXsession'};
        
        % Delete sessions that are all NaNs 
        cn_total =size(subject_patterns,3); 
        for cn = cn_total:-1:1 
            if nansum(nansum(subject_patterns(:,:,cn))) == 0
                subject_patterns(:,:,cn) = [];
            end
        end 
               
        % randomize trials, if needed
        if strcmp(p.Results.randomized_or_notrand, 'randomized') 
            subject_patterns = subject_patterns(:,:,randperm(size(subject_patterns,ndims(subject_patterns)))); %newly added
            allsubj_results.randomized = 'randomized';
        end
        
        % define the percentage of data that will be used as testing data
        if ~isfield(p.Results, 'test_percent') || isempty(p.Results.test_percent)
            test_percent = .2;
        else
            test_percent = p.Results.test_percent;
        end
           
        % define the folds
        fold_dim = ndims(subject_patterns);
        num_data = size(subject_patterns,fold_dim);
        num_in_fold = num_data*test_percent; 
        num_folds = floor(num_data/num_in_fold);
        num_in_fold = floor(num_in_fold); % round down 

    %% Leave-One-Out (loo)
    
    elseif strcmp(p.Results.approach, 'loo')   
        
        % take only one subject
        allsubj_results.approach = 'leave-one-out (loo)';
        
        if length(size(mcpa_summ.patterns)) == 5
            subject_patterns = mcpa_summ.patterns(:,:,:,:,s_idx);
        else
            subject_patterns = mcpa_summ.patterns(:,:,:,s_idx);
        end
        
        % randomize trials, if needed
        if strcmp(p.Results.randomized_or_notrand, 'randomized') 
            new_index = randperm(size(subject_patterns,ndims(subject_patterns)));
            subject_patterns = subject_patterns(:,:,new_index);
        end
        
        % define the fold
        num_folds = length(MCP_struct(s_idx).Experiment.Runs);
    end
    
    %% Begin cross-validation
    for folding_idx = 1:num_folds  
        
        %% Define fold_idx: indices in subject_patterns that will be test data
        temp_set_results_cond = nan(n_cond,n_sets,n_feature);
        
        if strcmp(p.Results.approach, 'kf')
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
        
        [train_data, train_labels, test_data, test_labels] = split_test_and_train(fold,...
            p.Results.conditions,...
            subject_patterns,... 
            mcpa_summ.event_types,...
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
        
        if strcmp(p.Results.approach, 'kf') && (~isfield(p.Results, 'n_randomsubset')||isempty(p.Results.n_randomsubset))
            allsubj_results.kf_n_randomsubset = 'no: using the whole dataset';
            allsubj_results.kf_sampling = 'no sampling';  
            
        elseif strcmp(p.Results.approach, 'kf') && isfield(p.Results, 'n_randomsubset')
            n_test = size(test_data,1);
            n_train = size(train_data,1);
            
            subset_test = floor(p.Results.n_randomsubset*test_percent);
            subset_train = floor(p.Results.n_randomsubset*(1-test_percent));
            allsubj_results.kf_n_randomsubset = sprintf('yes: %d out of %d test data and %d out of %d train data', subset_test, n_test, subset_train, n_train); 

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
        
        for set_idx = 1:min(n_sets,p.Results.max_sets)    

            if p.Results.verbose
                status_jump = floor(n_sets/20);
                if ~mod(set_idx,status_jump)
                    fprintf(' .')
                end
            end
            
            % Select the features for this subset
            set_features = sets(set_idx,:);

            %% Classify

            % RSA
            if strcmp(func2str(p.Results.test_handle),'rsa_classify')
                [predicted_labels, comparisons] = p.Results.test_handle(...
                    train_data(:,set_features,:), ...
                    train_labels,...
                    test_data(:,set_features,:),...
                    test_labels,...
                    p.Results.opts_struct);

            else
                [predicted_labels, comparisons] = p.Results.test_handle(...
                    train_data(:,set_features), ...
                    train_labels,...
                    test_data(:,set_features),...
                    test_labels,...
                    p.Results.opts_struct);
            end

            %% Record results 

            if size(predicted_labels,2) > 1 % test labels will be a column vector if we don't do pairwise
                if s_idx==1 && set_idx == 1 && folding_idx == 1
                    if strcmp(p.Results.approach, 'kf')   
                        allsubj_results.accuracy_matrix = nan(n_cond,n_cond,min(n_sets,p.Results.max_sets),num_folds,n_subj);                             
                    elseif strcmp(p.Results.approach, 'loo')   
                        allsubj_results.accuracy_matrix = nan(n_cond,n_cond,min(n_sets,p.Results.max_sets),n_sessions,n_subj); 
                    end
                end

                if iscell(comparisons)
                    subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
                    comparisons = cellfun(@(x) find(strcmp(x,mcpa_summ.event_types)),comparisons); 
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
                    test_labels(strcmp(strjoin(string(p.Results.conditions{cond_idx}),'+'),test_labels)),... % known labels
                    predicted_labels(strcmp(strjoin(string(p.Results.conditions{cond_idx}),'+'),test_labels))...% classifier labels
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
        
    if p.Results.verbose
        fprintf(' %0.1f mins\n',toc/60);
    end
    end % end session loop
    
end % end subject loop


%end % end function

