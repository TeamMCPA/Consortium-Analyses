function p = parse_inputs(MCP_struct, cv_function, varargin)
%% input parser 
% takes input from function that calls it and creates a struct to store
% clasification parameters
% if not input is provided for a particular variable, it creates default values

% Arguments:
% MCP_struct: MCP data structure
% varargin: input to main function

p = inputParser;
p.KeepUnmatched = true;

for s = 1:length(MCP_struct)
    % Some MCP files will not already have a transformation matrix
    % stored for translating channels into features. If that field is
    % missing, create an identity matrix. This field will be lost after
    % parsing the inputs, but the dimensions are needed to infer values for
    % incl_channels, incl_features, and setsize. If the transformation
    % matrix already exists, this default will not be used.
    if ~isfield(MCP_struct(s).Experiment.Runs(1),'Transformation_Matrix') || isempty(MCP_struct(s).Experiment.Runs(1).Transformation_Matrix)
        MCP_struct(s).Experiment.Runs(1).Transformation_Matrix = eye(length(MCP_struct(s).Experiment.Probe_arrays.Channels));
    end
end

%% parameters used for all kinds of classifiers
addParameter(p,'incl_channels',[1:max(arrayfun(@(x) size(x.Experiment.Runs(1).Transformation_Matrix,1),MCP_struct))],@isnumeric);
addParameter(p,'incl_subjects',[1:length(MCP_struct)],@isnumeric);
addParameter(p,'conditions',unique(cellstr(char(cellfun(@(x) char(x{:}), arrayfun(@(x) unique({x.Experiment.Conditions.Name},'stable'),MCP_struct, 'UniformOutput',false),'UniformOutput',false))),'stable'),@iscell);
addParameter(p,'summary_handle',@nanmean);
addParameter(p,'setsize',max(arrayfun(@(x) size(x.Experiment.Runs(1).Transformation_Matrix,2),MCP_struct)),@isnumeric);
addParameter(p,'max_sets',1000000,@isnumeric);
addParameter(p,'test_handle',@mcpa_classify);
addParameter(p,'opts_struct',[],@isstruct);
addParameter(p,'verbose',true);
addParameter(p, 'summarize_dimensions', {});
addParameter(p, 'final_dimensions', {});
addParameter(p, 'suppress_warnings', false);
addParameter(p, 'incl_features', [1:max(arrayfun(@(x) size(x.Experiment.Runs(1).Transformation_Matrix,2),MCP_struct))],@isnumeric); 


%% parameters that can have more than one value set
addParameter(p, 'hemoglobin', 'Oxy');
addParameter(p,'time_window',[2,6]);
addParameter(p,'baseline_window',[-5 0]);
addParameter(p,'scale_data', false);
addParameter(p,'incl_sessions',[1:max(arrayfun(@(x) length(x.Experiment.Runs),MCP_struct))]);

%% parameters specific to cross validation type
if strcmp(func2str(cv_function), 'classify_WithinSubjects')
    % parameters only used in nfold_classify_WithinSubjects
    addParameter(p, 'approach', 'loo', @ischar); 
    addParameter(p, 'randomized_or_notrand', 'notrand', @ischar); 
    addParameter(p, 'test_percent', []); % for kf 
    addParameter(p, 'n_randomsubset', [], @isnumeric); % for kf 
elseif strcmp(func2str(cv_function), 'generalize_ParticipantLevel')
    % parameters only used in nfold_generalize_ParticipantLevel
    addParameter(p,'cond_key',{});
    addParameter(p,'test_marks',{});
end

%% parameters used if norming the data
if any(cellfun(@(x) strcmp(x,'scale_data'), varargin(find(rem(1:length(varargin),2)))))
    addParameter(p,'scale_withinSessions', true, @islogical);
    addParameter(p, 'scale_function', @minMax_scale);
    addParameter(p, 'minMax', [0,1]);
end

%% parse the input and print out any unused parameters 
parse(p,varargin{:});


fields = fieldnames(p.Unmatched);
if ~isempty(fields)
    warning('The following fields were left unmatched: ')
    warning('%s ', fields{:})
end
  

end
