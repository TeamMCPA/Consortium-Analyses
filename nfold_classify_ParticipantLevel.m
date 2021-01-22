function allsubj_results = nfold_classify_ParticipantLevel(MCP_struct,varargin)
%% nfold_classify_ParticipantLevel takes an MCP struct and performs
% n-fold cross-validation for n subjects to classify individual
% participants' average response patterns.
%
%
% This wrapper allows the user to choose if
% features will be averaged within-participants to produce a single
% participant-level observation or if individual events will be preserved.
% If the featuers are averaged within-participants, the training set
% is constrained to the number of participants minus 1. Otherwise the
% training set will be (participants - 1) * (number of instances of an
% object).
%
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
% scale_data: option to perform some for of feature scaling. Defualt: false
% scale_withinSessions: if we do norm data, this allows the user to choose
% if we norm within individual sessions or across a participants session.
% It is only valid to norm across sessions for
% nfold_classify_ParticipantLevel. Default: true
% scale_function: function handle for how to scale the data. Default:
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
% Step 1: Epoching the data by time window

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
    isWithinSubjects = false;
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
    isWithinSubjects = false;
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

%% Begin the n-fold process: Select one test subj at a time from MCPA struct
for s_idx = 1:n_subj
    
    if input_struct.verbose
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

    [train_data, train_labels, test_data, test_labels] = split_test_and_train(s_idx,...
        input_struct.conditions,...
        mcpa_summ.patterns,...
        mcpa_summ.event_types,...
        final_dimensions,...
        mcpa_summ.dimensions, [], []);

    
    %% Run classifier and compare output with correct labels
    for set_idx = 1:min(n_sets,input_struct.max_sets)
        %% Progress reporting bit (not important to function. just sanity)
        % Report at every 5% progress
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
        %% Record results .

        if size(predicted_labels,2) > 1 % test labels will be a column vector if we don't do pairwise          
            subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
            nan_idx = cellfun(@(x) any(isnan(x)), predicted_labels(:,1,:), 'UniformOutput', false);
            subj_acc(:,:,[nan_idx{1,:,:}]) = nan;

            % Then loop through comparisons and save accuracy to the results struct
            for comp = 1:size(comparisons,1)
                allsubj_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),set_idx, s_idx) = subj_acc(comp);
            end

        else
            subj_acc = strcmp(predicted_labels, test_labels);
            for cond_idx = 1:n_cond
                cond_acc = nanmean(subj_acc(comparisons == cond_idx));
                allsubj_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = cond_acc;
                allsubj_results.accuracy(cond_idx).subjXfeature(s_idx,:) = cond_acc;
            end
        end
            

        %% Progress reporting
        if input_struct.verbose
            fprintf(' %0.1f mins\n',toc/60);
        end
    end % end set_idx loop
end % end subject loop

%% Visualization
if input_struct.verbose
    if n_sets > 1 && length(input_struct.conditions)==2
        
        figure
        errorbar(1:size(allsubj_results.accuracy(1).subjXfeature,2),mean(allsubj_results.accuracy(1).subjXfeature),std(allsubj_results.accuracy(1).subjXfeature)/sqrt(size(allsubj_results.accuracy(1).subjXfeature,1)),'r')
        hold;
        errorbar(1:size(allsubj_results.accuracy(2).subjXfeature,2),mean(allsubj_results.accuracy(2).subjXfeature),std(allsubj_results.accuracy(2).subjXfeature)/sqrt(size(allsubj_results.accuracy(2).subjXfeature,1)),'k')
        title('Decoding Accuracy across all features: Red = Cond1, Black = Cond2')
        set(gca,'XTick',[1:length(input_struct.incl_features)])
        set(gca,'XTickLabel',input_struct.incl_features)
        hold off;
        
        figure
        errorbar(1:size(allsubj_results.accuracy(1).subjXfeature,1),mean(allsubj_results.accuracy(1).subjXfeature'),repmat(std(mean(allsubj_results.accuracy(1).subjXfeature'))/sqrt(size(allsubj_results.accuracy(1).subjXfeature,2)),1,size(allsubj_results.accuracy(1).subjXfeature,1)),'r')
        hold;
        errorbar(1:size(allsubj_results.accuracy(2).subjXfeature,1),mean(allsubj_results.accuracy(2).subjXfeature'),repmat(std(mean(allsubj_results.accuracy(2).subjXfeature'))/sqrt(size(allsubj_results.accuracy(2).subjXfeature,2)),1,size(allsubj_results.accuracy(2).subjXfeature,1)),'k')
        title('Decoding Accuracy across all subjects: Red = Cond1, Black = Cond2')
        set(gca,'XTick',[1:input_struct.incl_subjects])
        hold off;
        
    end
end

end
