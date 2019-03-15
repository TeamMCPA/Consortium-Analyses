function [event_matrix] = MCP_get_subject_events(mcp_struct, channels, time_window, event_types)

%MCP_GET_SUBJECT_EVENTS Returns a matrix that contains HbO data for target 
%subject in each type, channel, time and type repetition.
%
% First we construct the index matrix that contains specific index of HbO
% data. Then Use the index matrix to find conresponding HbO data and add
% them into the output matrix.
%
% Chengyu Deng & Benjamin Zinszer 5 may 2017
% revised bdz 26 oct 2018

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
    max_condition_type = length(find(marks_vec(:, 1) == 1));
    
    for i = 2:length(event_types)
        if max_condition_type < length(find(marks_vec(:,i) == 1))
            max_condition_type = length(find(marks_vec(:, i) == 1));
        end
    end
    
    marks_mat = nan(max_condition_type, length(event_types));
    
    for type_i = 1:length(event_types)
        %Find the array index in the marks_vec
        temp_marks = find(marks_vec(:, type_i) == 1);

% This line was meant to remove every offset mark (every-other-mark) but it
% only applies if there *are* offsets, and if there aren't, it removes half
% of your events. Oops!
%         %Abandon the offsets
%         temp_marks = temp_marks(1:2:end);
        
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
time_window_samp = [round(time_window(1)*Fs_val) : round(time_window(end)*Fs_val)];
% The output matrix setup(time x channels x type repetition x types)
event_matrix = nan(length(time_window_samp), length(channels), size(marks_mat, 1), length(event_types));

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
        
        % For events where the full duration of time_window is available
        if marks_mat(event_j,type_i) + time_window_samp(1)>0 && ...
            marks_mat(event_j,type_i) + time_window_samp(end) <= length(oxy_timeser)
            
            % Grab the trial data (for this event) directly from the oxy
            % timeseries. It will be re-baselined momentarily.
            event_data = oxy_timeser(marks_mat(event_j,type_i)+time_window_samp,:);
            
            % Get the oxy data from one scan prior to stimulus onset and
            % use that as the baselining value.
            % There is good reason to think we should use maybe the average
            % of a few seconds prior to stimulus onset, but sometimes the 
            % events are only seconds apart, and we don't want to bleed
            % information about the previous trial into the current trial
            % by baselining off of it. For now, this approach is quick and
            % easy and probably wrong.
            %rebaseline_data = ones(length(time_window_samp),1)* ...
            %    oxy_timeser(marks_mat(event_j,type_i)-1,:);
            rebaseline_data = 0;

            event_matrix(:,:,event_j,type_i) = event_data - rebaseline_data;

            % This was the old way, and it caused problems. The new way is
            % easier to read, so you can see the problems more clearly.
            %event_matrix(:,:,event_j,type_i) = ...
            %    oxy_timeser(marks_mat(event_j,type_i)+time_window_samp,:) -...
            %    ones(length(time_window_samp),1)*oxy_timeser(marks_mat(event_j,type_i)+time_window_samp(1),:);
        
        % For events that might get truncated (e.g., end of recording)
        else
            % For the moment, we carelessly discard these trials. There is
            % a better way to do it, but you have to match up the lengths
            % of the vectors. Todo.
            event_matrix(:,:,event_j:end,type_i) = NaN;
        end
    end

end




end