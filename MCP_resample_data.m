

function new_mcp_file = MCP_resample_data(mcp_file,new_sampling_rate,save_flag)
%% Resample hemoglobin data as if the machine had a different sampling rate
% This is useful when doing timeXfeatures classification, or if comparing
% data with different sampling rates. 
% All resampling is done within session.
% Generates a new onsets, Oxy, Deoxy, and Total hemoglobin matrices.
% Updates indices for time sampling and for which samples belong to which
% session in the Runs field.

% input: mcp_file - can be the name of an MCP file or the MCP file itself
% new_sampling_rate - your desired sampling rate
% save_flag - save the new data locally


% Benjamin Zinszer & Anna Herbolzheimer 3 Sept 2020

%% Clean up / process inputs
% Open the old MCP struct or file
if isstruct(mcp_file)
    old_mcp_struct = mcp_file;
    if save_flag
        mcpfilename = ['MCP_dataset_' date];
        mcpdir = [pwd '/'];
        warning('File name not provided. Will save new data to working directory.')
    end
else
    [mcpdir, mcpfilename, ~] = fileparts(mcp_file);
    %fprintf('Loading file: %s\n\n',mcpfile);
    old_mcp_struct = load([mcpdir mcpfilename '.mcp'],'-mat');
end

if ~isnumeric(new_sampling_rate)
    error('new_sampling_rate must be a scalar numeric value, specified in Hz')
end

%% Recursion for multiple subjects (length MCP > 1)
% If there are multiple subjects, recurse the function for each one without
% saving and then save (if flagged) at the end
if length(old_mcp_struct) > 1
    for subj = 1:length(old_mcp_struct)
        disp(['Resampling data for: ' old_mcp_struct(subj).Subject.Subject_ID] );
        new_mcp_file(subj) = MCP_resample_data(old_mcp_struct(subj),new_sampling_rate,0);
    end
    if save_flag
        disp(['Saving new file under name: ' mcpfilename '_resampled.mcp']);
        disp([]);
        save([mcpdir mcpfilename '_resampled.mcp'],'-mat','new_mcp_file')
    end
    disp([]);
    return
end

%% Perform the resampling for one element of an MCP struct
new_mcp_file = old_mcp_struct;
new_mcp_file.fNIRS_Data.Hb_data.Oxy = [];
new_mcp_file.fNIRS_Data.Hb_data.Deoxy = [];
new_mcp_file.fNIRS_Data.Hb_data.Total = [];
new_mcp_file.fNIRS_Data.Onsets_Matrix = [];

% Obtain and check the old sampling rate
old_sampling_rate = new_mcp_file.fNIRS_Data.Sampling_frequency;
if old_sampling_rate < new_sampling_rate
    warning('The original data were sampled at a slower rate than the new sampling frequency. The data will be interpolated to increase sampling rate, but this is NOT recommended!');
end

if round(old_sampling_rate,4) ~= round(new_sampling_rate,4)
    
    % Each run is resampled separately to prevent interpolation between the
    % seaparate runs. Further, some runs may be sampled at different rates
    % than others.
    
    % Create a vector of old sampling rates (one per run) if that doesn't
    % already exist (based on the length of the old_sampling_rate).
    n_runs = length(new_mcp_file.Experiment.Runs);
    if length(old_sampling_rate)==1, old_sampling_rate = repmat(old_sampling_rate,n_runs,1); end
    
    % Step through each run and perform the resampling separately to
    % prevent data-leaking between runs
    
    idx_bookmark = 0; % Keep track of the new index values while concatenating data across runs
    
    for curr_run = 1:n_runs
        
        % Get the old data: frequency (sampling_rate), times, and indices
        old_Hz = old_sampling_rate(curr_run);
        old_Tx = old_mcp_struct.Experiment.Runs(curr_run).Time;
        old_Idx = old_mcp_struct.Experiment.Runs(curr_run).Index;
        
        % Use the indices for the run to extract the corresponding Hb data.
        old_oxy = old_mcp_struct.fNIRS_Data.Hb_data.Oxy(old_Idx,:);
        old_deoxy = old_mcp_struct.fNIRS_Data.Hb_data.Deoxy(old_Idx,:);
        old_total = old_mcp_struct.fNIRS_Data.Hb_data.Total(old_Idx,:);
        old_onset_matrix = old_mcp_struct.fNIRS_Data.Onsets_Matrix(old_Idx,:);

        % Generate the new Hb data
        new_Hz = new_sampling_rate;
        [new_oxy, new_Tx] = resample(old_oxy,old_Tx,new_Hz);
        [new_deoxy, ~] = resample(old_deoxy,old_Tx,new_Hz);
        [new_total, ~] = resample(old_total,old_Tx,new_Hz);
        
        % Generate the new onsets data.
        old_ons_idx = find(any(old_onset_matrix~=0,2));
        new_onsets = zeros(size(new_oxy,1),size(old_onset_matrix,2));
        % loop through the rows and find the corresponding time point in the
        % new onsets matrix
        for event_i = 1:length(old_ons_idx)
            old_event_Tx = old_Tx(old_ons_idx(event_i));
            [~, new_event_idx] = min(abs(new_Tx - old_event_Tx));
            new_onsets(new_event_idx,:) = old_onset_matrix(old_ons_idx(event_i),:);
        end
        
        % Generate the new run data
        new_Idx = [(idx_bookmark+1):1:(idx_bookmark+length(new_Tx))];
        new_mcp_file.Experiment.Runs(curr_run).Time = new_Tx;
        new_mcp_file.Experiment.Runs(curr_run).Index = new_Idx;
        
        % Update the idx_bookmark for the next run
        idx_bookmark = max(new_Idx);
        
        % Concatenate the Hb Data
        new_mcp_file.fNIRS_Data.Sampling_frequency(curr_run) = new_Hz;
        new_mcp_file.fNIRS_Data.Hb_data.Oxy(new_Idx,:) = new_oxy;
        new_mcp_file.fNIRS_Data.Hb_data.Deoxy(new_Idx,:) = new_deoxy;
        new_mcp_file.fNIRS_Data.Hb_data.Total(new_Idx,:) = new_total;
        new_mcp_file.fNIRS_Data.Onsets_Matrix(new_Idx,:) = new_onsets;
        
    end

else
    warning('Old and new sampling rates equal to the 0.0001 decimal place. Taking no action.')
    new_mcp_file = old_mcp_struct;
end

% If all runs are resampled to same frequency (sampling rate; and they
% should be by now!) then reduce the vector to a scalar.
if length(unique(new_mcp_file.fNIRS_Data.Sampling_frequency))==1
    new_mcp_file.fNIRS_Data.Sampling_frequency = unique(new_mcp_file.fNIRS_Data.Sampling_frequency);
end

%% If the save_flag is true, write the data out.
if save_flag, save([mcpdir mcpfilename '_resampled.mcp'],'-struct','new_mcp_file'); end


