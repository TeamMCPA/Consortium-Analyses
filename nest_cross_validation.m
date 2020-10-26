function final_results = nest_cross_validation(MCP_struct, cv_function, varargin)
%% Load MCP struct if necessary
if isstring(MCP_struct) || ischar(MCP_struct)
    MCP_struct = load(MCP_struct,'-mat');
    varname = fieldnames(MCP_struct);
    MCP_struct = eval(['MCP_struct.' varname{1}]);
end

%% Parse out the input data
p = parse_inputs(MCP_struct, cv_function, varargin{:});   

%% Setting up the combinations of feature subsets
unmapped_sets = find_feature_sets(p.Results);
sets = map_features_to_sets(p, unmapped_sets);

%% establish how to concatenate or average over our dimensions
% intermediary step: see if the user specified the summarizing dimensions. If not,
% recommend what dimensions to average over

if ~isempty(p.Results.summarize_dimensions) || ~isfield(p.Results, 'summarize_dimensions')
    summarize_dimensions = p.Results.summarize_dimensions;
else
    isWithinSubjects = false;
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
    isWithinSubjects = false;
    warning('final_dimensions not specified. Consulting recommend_dimensions.')
    
    [~, final_dimensions] = recommend_dimensions(p.Results, isWithinSubjects);

    fprintf('The format the data will be in when it enters the classifier wrapper is: %s', final_dimensions{:});
    fprintf('\n')
end

%% Prep basic parameters
n_subj = length(p.Results.incl_subjects);
n_sets = size(sets,1);
n_feature = length(p.Results.incl_features);
try n_cond = length(unique(p.Results.conditions)); catch, n_cond = length(p.Results.conditions); end

%% create results struct
final_results = create_results_struct(p.Results,...
    cv_function,...
    sets,...
    length(p.Results.incl_subjects),...
    size(sets,1),...
    length(p.Results.incl_features),...
    length(p.Results.conditions),...
    [],...
    final_dimensions,...
    [],...
    []);

[parameter_space, proc_params, classif_params] = enumerate_parameter_space(final_results);

%% begin the cross validation
for s_idx = 1:n_subj
    %% test the different parameters

    train_data_idx = setdiff(1:n_subj, s_idx);
    nesting_results = struct;

    for p_idx = 1:size(parameter_space,1)
        inner_results_struct = populate_classification_options_struct(final_results,... 
            parameter_space,...
            proc_params,...
            classif_params,...
            p_idx);
        
        inner_results_struct.incl_subjects = 1:(length(inner_results_struct.incl_subjects)-1);
        
        if inner_results_struct.verbose
            inner_results_struct.verbose = false;
            fprintf('running nested cross validation for parameter set: %g / %g', p_idx, size(parameter_space,1))
            fprintf('\n')
        end
        
        % create mcpa struct
        mcpa_struct = MCP_to_MCPA(MCP_struct(train_data_idx),...
            (1:length(train_data_idx)),...
            inner_results_struct.incl_features,...
            inner_results_struct.incl_channels,...
            inner_results_struct.time_window,...
            inner_results_struct.baseline_window,...
            inner_results_struct.hemoglobin);

        % Subset patterns by session
        if max(inner_results_struct.incl_sessions) > size(mcpa_struct.patterns,ndims(mcpa_struct.patterns)-1)
            inner_results_struct.incl_sessions = size(mcpa_struct.patterns,ndims(mcpa_struct.patterns)-1);
        end
        inds = pad_dimensions(mcpa_struct.dimensions, 'session', inner_results_struct.incl_sessions);
        mcpa_struct.patterns = mcpa_struct.patterns(inds{:}); % Replace patterns matrix with the subsetted sessions data

        % summarize the data
        mcpa_summ = summarize_MCPA_Struct(inner_results_struct.summary_handle,...
            mcpa_struct,...
            summarize_dimensions);
        
        inner_results_struct.patterns = mcpa_summ.patterns;
        inner_results_struct.event_types = mcpa_summ.event_types;
        inner_results_struct.dimensions = mcpa_summ.dimensions;
        
        inner_results_struct.summed_mcpa_patterns = mcpa_summ.patterns;
        inner_results_struct.summarize_dimensions = mcpa_summ.summarize_dimensions;
        
        inner_results_struct = cv_function(inner_results_struct); 
        
        nesting_results.(['Parameter_set_' int2str(p_idx)]).accuracy = inner_results_struct.accuracy;
        nesting_results.(['Parameter_set_' int2str(p_idx)]).parameters = parameter_space{p_idx,:};
        if isfield(inner_results_struct, 'accuracy_matrix')
            nesting_results.(['Parameter_set_', int2str(p_idx)]).accuracy_matrix = inner_results_struct.accuracy_matrix;
        end

        
    end
       
    %% evalute the results
    [chosen_parameters, rankings] = eval_max_accuracy(nesting_results, parameter_space);
    nesting_results.('chosen_parameters') = chosen_parameters;
    nesting_results.('rankings') = rankings;
    
    %% apply to left out data
    outer_results_struct = populate_classification_options_struct(final_results,...
        parameter_space,...
        proc_params,...
        classif_params,...
        rankings(1));

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
        summarize_dimensions);
    outer_results_struct.event_types = mcpa_summ.event_types;
    
    % split the data
    [train_data, train_labels, test_data, test_labels] = split_test_and_train(s_idx,...
        p.Results.conditions,...
        mcpa_summ.patterns,...
        mcpa_summ.event_types,...
        final_dimensions,...
        mcpa_summ.dimensions, [], []);

    for set_idx = 1:min(n_sets,outer_results_struct.max_sets)
        %% Progress reporting bit (not important to function. just sanity)
        % Report at every 5% progress
        if p.Results.verbose
            status_jump = floor(n_sets/20);
            if ~mod(set_idx,status_jump)
                fprintf(' .')
            end
        end
        % Select the features for this subset
        set_features = sets(set_idx,:);
        

        %% Classify
        inds = pad_dimensions(final_dimensions, 'feature', set_features);
        [predicted_labels, comparisons] = outer_results_struct.test_handle(...
                train_data(inds{:}), ...
                train_labels,...
                test_data(inds{:}),...
                test_labels,...
                outer_results_struct.opts_struct);
            
           
        outer_results_struct = save_accuracy(outer_results_struct,...
            comparisons,...
            predicted_labels,...
            set_idx,...
            1,... 
            s_idx);    
            
        
        
    end
    
    nesting_results.final_accuracy = outer_results_struct.accuracy;
    if isfield(outer_results_struct, 'accuracy_matrix')
        nesting_results.final_accuracy_matrix = outer_results_struct.accuracy_matrix;
    end
    
    final_results.Nesting_results(s_idx).Results = nesting_results;

    
end


end