function [event_matrix] = MCP_get_subject_events(mcp_struct, channels, time_window, event_types, base_window)

%MCP_GET_SUBJECT_EVENTS Returns a matrix that contains HbO data for target 
%subject in each type, channel, time and type repetition.
%
% First we construct the index matrix that contains specific index of HbO
% data. Then Use the index matrix to find conresponding HbO data and add
% them into the output matrix.
%
% Chengyu Deng & Benjamin Zinszer 5 may 2017
% revised bdz 29 aug 2019

%% If no base-lining window (base_window) is specified, mark NaN flag
if ~exist('base_window','var')
    base_window = NaN;
end

%% MCP_get_subject_events will run recursively for multiple subjects
% If the mcp_struct contains multiple elements (i.e., a struct for each
% subject), then the function will call itself and re-evaluate on the
% contents of each struct, returning event_matrix in a separate cell for
% each element in the struct (i.e., a cell for each subject).
if length(mcp_struct)>1
    event_matrix = cell(length(mcp_struct),1);
    for subj_num = 1:length(mcp_struct)
        event_matrix{subj_num} = MCP_get_subject_events(mcp_struct(subj_num), channels, time_window, event_types);
    end
    return
end

%% Index matrix handling (We use MCP struct so we don't need distinguish Homer file version anymore)

% Extract oxy data and marks from the MCP struct
oxy_timeser = mcp_struct.fNIRS_Data.Hb_data.Oxy(:, channels);
marks_vec = mcp_struct.fNIRS_Data.Onsets_Matrix;

% Handle different type of marks vector

if size(marks_vec, 2) > 1
    
    % Determine the maximum number of reps for any given marker (iterates
    % through the columns of marks_vec and finds values >0, aka events)
    max_condition_type = max(sum(marks_vec>0,1));
%     max_condition_type = 0;
%     for i = 1:size(marks_vec,2)
%         if max_condition_type < length(find(marks_vec(:,i) == 1))
%             max_condition_type = length(find(marks_vec(:, i) == 1));
%         end
%     end
    
    marks_mat = nan(max_condition_type, size(marks_vec,2));
    
    for type_i = 1:size(marks_vec,2)
        %Find the array index in the marks_vec
        temp_marks = find(marks_vec(:, type_i) == 1); 
        marks_mat(1:length(temp_marks), type_i) = temp_marks;
    end
    
    % Remove rows whose entries are all NaN
    marks_mat = marks_mat(sum(isnan(marks_mat),2)<size(marks_mat,2),:);
    
elseif size(marks_vec, 2) == 1
    
    % marks_mat matrix's row is the max number of the condition in whole
    % conditions, and the column is the event types.
    marks_mat = nan(max(hist(marks_vec(marks_vec~=0))),length(event_types));
    
    for type_i = 1 : length(event_types)
        %Find the array index (vector index) in the marks_vec
        temp_marks = find(marks_vec == event_types(type_i));
        
% This line was meant to remove every offset mark (every-other-mark) but it
% only applies if there *are* offsets, and if there aren't, it removes half
% of your events. Oops!
        %Now abandon the offsets
%        temp_marks = temp_marks(1:2:end);
        
        marks_mat(1:length(temp_marks), type_i) = temp_marks;
    end
    
    % Remove rows whose entries are all NaN
    marks_mat = marks_mat(sum(isnan(marks_mat),2)<size(marks_mat,2),:);
    
else
    %We shall never reach this part hopefully.
    fprintf('Unidentified data format. Fail to return right event matrix.');
    event_matrix = [];
    return
end

%% Extract the individual events for each event type
Fs_val = mcp_struct.fNIRS_Data.Sampling_frequency;
time_window_samp = round(time_window.*Fs_val);
base_window_samp = round(base_window.*Fs_val);
% The output matrix setup(time x channels x type repetition x types)
num_samps = max(time_window_samp) - min(time_window_samp) + 1;
event_matrix = nan(num_samps, length(channels), size(marks_mat, 1), length(event_types));

for type_i = 1 : length(event_types)
    
    matched_conditions = cell2mat(arrayfun(@(x) strcmp(event_types{type_i},x.Name), mcp_struct.Experiment.Conditions,'UniformOutput',false)); 
    event_marks = [mcp_struct.Experiment.Conditions(matched_conditions).Mark];
    
    % Each instance of a given condition (i.e., each individual event) is
    % extracted and entered into the event_matrix. Since there are often an
    % unbalanced number of events per condition, the events_matrix is
    % preallocated to accomodate the largest number of events for any given
    % condition, and the conditions with fewer events have NaN values
    % filling in the rest. These NaNs are important to remember if you're
    % taking a mean value later (use mean(x,dim,'omitnan) or nanmean(x,dim)
    for event_j = 1 : length(marks_mat(:, event_marks))
        
        % Find the earliest time-point (in samples) that will be extracted
        % and the latest time-point (in samples) that will be extracted so
        % that we can determine the window of data to be extracted.
        earliest_samp = min( [time_window_samp(:); base_window_samp(:)] );
        latest_samp = max( [time_window_samp(:); base_window_samp(:)] );
        
        % For events where the full duration of the baseline+event is 
        % available, extract the event_data (within the time_window) and
        % rebaseline_data (within the baseline_window). rebaseline_data are
        % then averaged, and the average is removed from the event_data
        if marks_mat(event_j,type_i) + earliest_samp > 0 && ...
            marks_mat(event_j,type_i) + latest_samp <= length(oxy_timeser)
            
            % Grab the trial data (for this event) directly from the oxy
            % timeseries. It will be re-baselined momentarily.
            event_data = oxy_timeser(marks_mat(event_j,type_i)+min(time_window_samp):marks_mat(event_j,type_i)+max(time_window_samp),:);
            
            % Get the oxy data from one scan prior to stimulus onset and
            % use that as the baselining value.
            if ~isnan(base_window)
                rebaseline_data = oxy_timeser(marks_mat(event_j,type_i)+min(base_window_samp):marks_mat(event_j,type_i)+max(base_window_samp),:);
            else
                rebaseline_data = 0;
            end
            event_matrix(:,:,event_j,type_i) = event_data - nanmean(rebaseline_data);

            % This was the old way, and it caused problems. The new way is
            % easier to read, so you can see the problems more clearly.
            %event_matrix(:,:,event_j,type_i) = ...
            %    oxy_timeser(marks_mat(event_j,type_i)+time_window_samp,:) -...
            %    ones(length(time_window_samp),1)*oxy_timeser(marks_mat(event_j,type_i)+time_window_samp(1),:);
        
        % For events that might get truncated (start or end of recording)
        else
            % For the moment, we carelessly discard these trials. There is
            % probably a better way to do it, but you have to match up the 
            % lengths of the vectors and figure out rebaselining. Todo.
            event_matrix(:,:,event_j:end,type_i) = NaN;
        end
    end

end




end