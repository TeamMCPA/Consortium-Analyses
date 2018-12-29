function MCPA_struct = MCP_to_MCPA(mcp_multiple, incl_subjects, incl_channels, time_window)
%MCP_TO_MCPA Convert MCP format data to MCPA_struct for analysis
% The function is called with the following arguments:
% MCP_to_MCPA(mcp_struct, incl_subjects, incl_channels, time_window)
%
% mcp: An customized MCP struct that contains all data for the analysis.
% Using MCP struct to store data can unify the way that data stored in the
% struct. Directly grabbing data from homer file might cause problems such
% as failure to find specific data.
%
% incl_subjects: a vector of indices for subjects to include in the
% analysis. Importantly the subject numbers correspond to the index in the
% struct array (e.g., MyData([1 3 5]) not any other subject number
% assignment. Use [] to just get all subjects.
%
% incl_channels: a vector of indices for channels to include in the
% analysis. Again, only the channel's position in the HomER struct matters,
% not any other channel number assignment. Use [] to get all channels.
%
% time_window: defined in number of seconds. If two subjects have different
% sampling frequencies, the same time window will be searched (except for
% rounding error of first and last samples). Time window can be specified
% as either [start, end] or [start : end] since only first and last times
% are used. Default (use []) is [-5,20] sec.
%
% The function will return a new struct containing some metadata and the
% multichannel patterns for each participant and condition.
%
% Chengyu Deng & Benjamin Zinszer 5 may 2017
% revised bdz 26 oct 2018

%% Check whether importing an MCP file or just converting from workspace
% Pulling from a file will be much faster for individual event
% classification later, so this method is preferred.

if isstruct(mcp_multiple)
    no_mcp_file = true;
else
    if iscell(mcp_multiple), mcp_multiple = mcp_multiple{:}; end
    MCPA_struct.data_file = {mcp_multiple};
    no_mcp_file = false;
    mcp_file_content = load(mcp_multiple,'-mat');
    varname = fieldnames(mcp_file_content);
    mcp_multiple = eval(['mcp_file_content.' varname{1}]);
    clear('mcp_file_content')
end

%% Double-check for missing data
if ~exist('incl_subjects','var') || isempty(incl_subjects)
    incl_subjects = 1:length(mcp_multiple);
end
if ~exist('incl_channels','var') || isempty(incl_channels)
    incl_channels = [1:max(arrayfun(@(x) size(x.fNIRS_Data.Hb_data.Oxy,2),mcp_multiple))];
end
if ~exist('time_window','var') || isempty(time_window)
    time_window = [-5,20];
end

%% Convert time window from seconds to scans
% rounds off sampling frequencies to 8 places to accomodate floating point
% errors
Fs_val = unique(round(arrayfun(@(x) x.fNIRS_Data.Sampling_frequency,mcp_multiple),8));
if length(Fs_val) > 1
    minFs = min(Fs_val);
    maxFs = max(Fs_val);
    warning([num2str(length(Fs_val)) ' different sampling frequencies found, ranging from ' num2str(minFs) ' to ' num2str(maxFs) ' Hz. '...
        'Data may be resampled during analysis if comparisons are made in time-domain']);
end
num_time_samps = length(round(time_window(1)*max(Fs_val)) : round(time_window(end)*max(Fs_val)));


%% Event type Handling

% % This version just gets integers up to the max number of conditions
% event_types = 1:max(arrayfun(@(x) length(x.Experiment.Conditions),MCP_data));

% This version uses names from the MCP array
all_names = arrayfun(@(x) unique({x.Experiment.Conditions.Name}),mcp_multiple(incl_subjects), 'UniformOutput',false);
unique_names = cellfun(@(x) char(x{:}),all_names,'UniformOutput',false);
event_types = unique(cellstr(char(unique_names)));

%% Extract data from the data file into the empty output matrix

% Initiate the subj_mat matrix that will be output later(begin with NaN)
subj_mat = nan(num_time_samps, length(event_types), length(incl_channels), length(incl_subjects));
fprintf('Output matrix for MCPA_struct is in dimension: time_window x types x channels x subjects\n');

% Extract data from each subject
fprintf('\nExtracting data for subject: \n');

for subj_idx = 1 : length(incl_subjects)
    
    fprintf('subject: %d. \n', incl_subjects(subj_idx));
    
    if no_mcp_file
    MCPA_struct.data_file{subj_idx} = [mcp_multiple(incl_subjects(subj_idx)).Experiment.Runs.Source_files]';
    end
    
    % Event_matrix format:
    % (time x channels x repetition x types)
    event_matrix = MCP_get_subject_events(mcp_multiple(incl_subjects(subj_idx)), incl_channels, time_window, event_types);
    
    % Event_repetition_mean:
    % (time x channels x event_type)
    event_repetition_mean = nanmean(event_matrix, 3);
    event_repetition_mean = reshape(event_repetition_mean, size(event_matrix,1), size(event_matrix,2), size(event_matrix,4));
    event_repetition_mean = permute(event_repetition_mean, [1 3 2]);
    % Now the dimension of Event_repetition_mean:
    % (time x event_types x channels )
    
    subj_time_samps = size(event_repetition_mean,1);
    if subj_time_samps == num_time_samps,
        % Output format: subj_mat(time_window x event_types x channels x subjects)
        subj_mat(:, :, :, subj_idx) = event_repetition_mean;
    else
        time_mask = round(linspace(1,num_time_samps,subj_time_samps));
        % Output format: subj_mat(time_window x event_types x channels x subjects)
        subj_mat(time_mask, :, :, subj_idx) = event_repetition_mean;
        
    end
    
end

%% Return the MCPA_struct
fprintf('\nWriting MCPA_struct for this dataset...');
try
    MCPA_struct.created = datestr(now);
    MCPA_struct.time_window = [round(time_window(1)*max(Fs_val))/max(Fs_val) : 1/max(Fs_val) : round(time_window(end)*max(Fs_val))/max(Fs_val)];
    MCPA_struct.incl_subjects = incl_subjects;
    MCPA_struct.incl_channels = incl_channels;
    MCPA_struct.event_types = event_types;
    MCPA_struct.patterns = subj_mat;
    
    fprintf('Done.\n');
catch
    fprintf('Failed to create the new struct (MCP_to_MCPA).\n');
    MCPA_struct = struct;
end


end