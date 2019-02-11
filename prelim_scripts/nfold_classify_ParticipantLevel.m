function allsubj_results = nfold_classify_ParticipantLevel(MCP_struct,varargin)
%% nfold_classify_ParticipantLevel takes an MCP struct and performs
% n-fold cross-validation for n subjects to classify individual 
% participants. This wrapper assumes that features will be
% averaged within-subjects to produce a single subject-level observation.
% Thus the training set is constrained to the number of subjects minus 1.
% Several parameters can be changed, including which functions are used to
% generate features and what classifier is trained. See Arguments below:
%
% Arguments:
% incl_channels: channels to include in the analysis. Default: all channels
% incl_subjects: index of subjects to include. Default: all subjects
% time_window: [onset, offset] in seconds. Default [2,6]
% conditions: cell array of condition names / trigger #s. Default: {1,2}
% summary_handle: function handle (or char of function name) to specify how
% time-x-channel data should be summarized into features. Default: nanmean
% setsize: number of channels to analyze (for subset analyses) Default: all
% test_handle: function handle for classifier. Default: mcpa_classify
% opts_struct: contains additional classifier options. Default: empty struct

%% Load MCP struct if necessary
if isstring(MCP_struct)
    MCP_struct = load(MCP_struct,'-mat');
    varname = fieldnames(MCP_struct);
    MCP_struct = eval(['MCP_struct.' varname{1}]);
end

%% Parse out the input data
p = inputParser;
addParameter(p,'incl_channels',[1:max(arrayfun(@(x) size(x.fNIRS_Data.Hb_data.Oxy,2),MCP_struct))],@isnumeric);
addParameter(p,'incl_subjects',[1:length(MCP_struct)],@isnumeric);
addParameter(p,'time_window',[2,6],@isnumeric);
addParameter(p,'conditions',{1,2},@iscell);
addParameter(p,'summary_handle',@nanmean);
addParameter(p,'setsize',max(arrayfun(@(x) size(x.fNIRS_Data.Hb_data.Oxy,2),MCP_struct)),@isnumeric);
addParameter(p,'test_handle',@mcpa_classify);
addParameter(p,'opts_struct',struct,@isstruct);
parse(p,varargin{:})

%% Build MCPA struct for all subjects in the MCP

% Step 1: Epoching the data by time window and averaging the epochs
% together at the subject level
mcpa_struct = MCP_to_MCPA(MCP_struct,p.Results.incl_subjects,p.Results.incl_channels,p.Results.time_window);

% Step 2: Apply the desired function (e.g., @nanmean) for summarizing time
% window data. You can write custom functions to deal with time- and
% channel-domain data however you want. Default behavior is to apply the
% function along the first dimension of the MCPA pattern, but this can also
% be changed.
mcpa_summ = summarize_MCPA_Struct(p.Results.summary_handle,mcpa_struct);


if length(p.Results.conditions)==2
    
    % This is the special case where there are only two conditions to
    % discriminate, and RSA cannot be applied. Instead, a two-class
    % classifier is applied directly (such as MCPA or SVM) to the data.
    cond1 = p.Results.conditions{1};
    cond2 = p.Results.conditions{2};
    
    %% Identify the conditions that will be compared
    % Condition flags can be string, cell, integer, or logical. If string or
    % cell, a search is run over the summarized MCPA data to determine which
    % condition # matches the condition name. Otherwise, integer or logical
    % values are applied directly to indexing.
    if ischar(cond1) || isstring(cond1) || iscellstr(cond1)
        cond1_flag = strcmp(cond1,mcpa_summ.event_types);
    else
        cond1_flag = cond1;
    end
    if ischar(cond2) || isstring(cond2) || iscellstr(cond2)
        cond2_flag = strcmp(cond2,mcpa_summ.event_types);
    else
        cond2_flag = cond2;
    end
    
    allsubj_results = leave_one_Ss_out_classifyAverages(mcpa_summ,cond1_flag,cond2_flag,p.Results.setsize);
    
else
    
    % Write the multiclass version here
    allsubj_results = pairwise_rsa_leaveoneout(mcpa_summ.patterns);
    
end

end