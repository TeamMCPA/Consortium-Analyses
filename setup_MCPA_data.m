function allsubj_results = setup_MCPA_data(MCP_struct,varargin)
%% Parse out the input data
input_struct = parse_inputs(MCP_struct, varargin{1}{:});  

%% assess what kind of cross validation this is
stack = dbstack;
current_folding_function = stack(2).name;

%% validate classification options
input_struct.opts_struct = validate_classification_options_input(MCP_struct, input_struct, input_struct.suppress_warnings);

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
    isWithinSubjects = strcmp(current_folding_function, 'nfold_classify_WithinSubjects');
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
    isWithinSubjects = strcmp(current_folding_function, 'nfold_classify_WithinSubjects');
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
allsubj_results = create_results_struct(input_struct,...
    current_folding_function,...
    sets,...
    n_subj,...
    n_sets,...
    n_feature,...
    n_cond,...
    mcpa_summ.dimensions,...
    final_dimensions,...
    mcpa_summ.patterns,...
    mcpa_summ.event_types);
    



end