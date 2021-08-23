function final_results = nest_nfold_classify_ParticipantLevel(MCP_struct, varargin)
%% performs a nested cross validation for between subjects comparisons
% nest_nfold_classify_ParticipantLevel optimizes parameters for an n-fold cross validation,
% where on each fold, one subject is held out as test data, and the rest
% are used to create a group model. On each fold, before classifiying the
% test data, a second, inner cross validation is performed on the training data only, where all combinations
% of user specfified parameters are tested to find the one combination that leads to the
% highest accuracy. That one combination of parameters will then be used to
% classify the held out test subject in the outer cross validation.

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
input_struct= parse_inputs(MCP_struct, varargin{:});   

%% Setting uinput_structthe combinations of feature subsets
unmapped_sets = find_feature_sets(input_struct);
sets = map_features_to_sets(input_struct, unmapped_sets);

%% establish how to concatenate or average over our dimensions
% intermediary step: see if the user specified the summarizing dimensions. If not,
% recommend what dimensions to average over

if ~isempty(input_struct.summarize_dimensions) || ~isfield(input_struct, 'summarize_dimensions')
    summarize_dimensions = input_struct.summarize_dimensions;
else
    warning('summarize_dimensions not specified. Consulting recommend_dimensions.')
    isWithinSubjects = false;
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


%% Prepare basic parameters
n_subj = length(input_struct.incl_subjects);
n_sets = size(sets,1);
n_feature = length(input_struct.incl_features);
try n_cond = length(unique(input_struct.conditions)); catch, n_cond = length(input_struct.conditions); end

%% create results struct
final_results = create_results_struct(input_struct,...
    'nfold_classify_ParticipantLevel',...
    sets,...
    length(input_struct.incl_subjects),...
    size(sets,1),...
    length(input_struct.incl_features),...
    length(input_struct.conditions),...
    summarize_dimensions,...
    final_dimensions,...
    [],...
    []);


%% set up a cell array with all parameter combinations 
[parameter_space, proc_params, classif_params] = enumerate_parameter_space(final_results);

%% begin the cross validation
for s_idx = 1:n_subj 
    % initialize a storage for inner optimization results
    nesting_results=struct; 
    inner_results = cell(size(parameter_space,1), 1);
    inner_options_structs = cell(size(parameter_space,1), 1);
    
    % set up data for inner optimization
    train_data_idx = setdiff(1:n_subj, s_idx); 
    
    % progress reporting
    if final_results.verbose == 3
        fprintf('running nested cross validation for subject %g / % g', s_idx, n_subj)
        fprintf('\n')
    end

    % begin optimization
    for p_idx = 1:size(parameter_space,1)
       
        % populate a results struct containing this optimization loop's
        % parameters
        inner_results_struct = populate_classification_options_struct(final_results,... 
            parameter_space,...
            proc_params,...
            classif_params,...
            p_idx,...
            true);

        % progress reporting
        if inner_results_struct.verbose == 1 || inner_results_struct.verbose == 2
            fprintf('running nested cross validation for subject %g / % g and parameter set: %g / %g', s_idx, n_subj, p_idx, size(parameter_space,1))
            fprintf('\n')
        end

        % optional - scale the data
        if inner_results_struct.scale_data
            MCP_struct2 = scale_individuals(MCP_struct, inner_results_struct);
        else
            MCP_struct2 = MCP_struct;
        end

        % create mcpa struct
        mcpa_struct = MCP_to_MCPA(MCP_struct2(train_data_idx),...
            (1:length(train_data_idx)),...
            inner_results_struct.incl_features,...
            inner_results_struct.incl_channels,...
            inner_results_struct.time_window,...
            inner_results_struct.baseline_window,...
            inner_results_struct.hemoglobin);

        % Subset patterns by session
        if max(inner_results_struct.incl_sessions) > size(mcpa_struct.patterns,ndims(mcpa_struct.patterns)-1)
            inner_results_struct.incl_sessions = 1:size(mcpa_struct.patterns,ndims(mcpa_struct.patterns)-1);
        end
        inds = pad_dimensions(mcpa_struct.dimensions, 'session', inner_results_struct.incl_sessions);
        mcpa_struct.patterns = mcpa_struct.patterns(inds{:}); % Replace patterns matrix with the subsetted sessions data

        % summarize the data
        mcpa_summ = summarize_MCPA_Struct(inner_results_struct.summary_handle,...
            mcpa_struct,...
            summarize_dimensions);

        % update inner results struct with newly created parameters for
        % classification
        inner_results_struct.patterns = mcpa_summ.patterns;
        inner_results_struct.event_types = mcpa_summ.event_types;
        inner_results_struct.dimensions = mcpa_summ.dimensions;
        inner_results_struct.summarize_dimensions = mcpa_summ.summarize_dimensions;

        % run the inner loop
        inner_results_struct = nfold_classify_ParticipantLevel(MCP_struct2,...
            'results_struct', inner_results_struct); 
        
        % store results
        nesting_results.(['Parameter_set_' int2str(p_idx)]).accuracy = inner_results_struct.accuracy;
        nesting_results.(['Parameter_set_' int2str(p_idx)]).parameters = parameter_space{p_idx,:};
        if isfield(inner_results_struct, 'accuracy_matrix')
            nesting_results.(['Parameter_set_', int2str(p_idx)]).accuracy_matrix = inner_results_struct.accuracy_matrix;
        end

    end
    
    %% test how well the optimization did on the test data
    %% evaluate which parameters were best
    [chosen_parameters, rankings] = eval_max_accuracy(nesting_results, parameter_space);
    nesting_results.('chosen_parameters') = chosen_parameters;
    nesting_results.('rankings') = rankings;

    %% apply to left out data
    outer_results_struct = populate_classification_options_struct(final_results,...
        parameter_space,...
        proc_params,...
        classif_params,...
        rankings(1),...
        false);
    
    % validate classification options
    outer_results_struct.opts_struct = validate_classification_options_input(MCP_struct, outer_results_struct, outer_results_struct.suppress_warnings);
    
    % scale the data
    if outer_results_struct.scale_data
        MCP_struct = scale_individuals(MCP_struct, outer_results_struct);
    end

    
    % create mcpa struct
    mcpa_struct = MCP_to_MCPA(MCP_struct,...
        outer_results_struct.incl_subjects,...
        outer_results_struct.incl_features,...
        outer_results_struct.incl_channels,...
        outer_results_struct.time_window,...
        outer_results_struct.baseline_window,...
        outer_results_struct.hemoglobin);

    % Subset patterns by session
    inds = pad_dimensions(mcpa_struct.dimensions, 'session', outer_results_struct.incl_sessions);
    mcpa_struct.patterns = mcpa_struct.patterns(inds{:}); % Replace patterns matrix with the subsetted sessions data

    % summarize the data
    mcpa_summ = summarize_MCPA_Struct(outer_results_struct.summary_handle,...
        mcpa_struct,...
        inner_results_struct.summarize_dimensions);
    outer_results_struct.event_types = mcpa_summ.event_types;
    
    % split the data
    [train_data, train_labels, test_data, test_labels] = split_test_and_train(s_idx,...
        input_struct.conditions,...
        mcpa_summ.patterns,...
        mcpa_summ.event_types,...
        inner_results_struct.final_dimensions,...
        mcpa_summ.dimensions, [], []);

    for set_idx = 1:min(n_sets,outer_results_struct.max_sets)
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
        inds = pad_dimensions(inner_results_struct.final_dimensions, 'feature', set_features);
        [predicted_labels, comparisons] = outer_results_struct.test_handle(...
                train_data(inds{:}), ...
                train_labels,...
                test_data(inds{:}),...
                test_labels,...
                outer_results_struct.opts_struct);
                            
        %% Record results 
        if size(predicted_labels,2) > 1 % test labels will be a column vector if we don't do pairwise          
            subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
            nan_idx = cellfun(@(x) any(isnan(x)), predicted_labels(:,1,:), 'UniformOutput', false);
            subj_acc(:,:,[nan_idx{1,:,:}]) = nan;

            % Then loop through comparisons and save accuracy to the results struct
            for comp = 1:size(comparisons,1)
                final_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),set_idx, s_idx) = subj_acc(comp);
            end

        else
            subj_acc = double(strcmp(predicted_labels, test_labels));
            nan_idx = cellfun(@isnan, predicted_labels);
            subj_acc(nan_idx) = NaN;
            for cond_idx = 1:n_cond
                cond_acc = nanmean(subj_acc(strcmp(comparisons, test_labels(cond_idx))));
                final_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = cond_acc;
                final_results.accuracy(cond_idx).subjXfeature(s_idx,:) = cond_acc;
            end
        end        
        
    end
    
    %% update the results struct with nesting results
    
    final_results.Nesting_results(s_idx).Optimization_Results = nesting_results;

    
end


end


