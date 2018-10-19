function [event_matrix] = MCP_get_subject_events(mcp_multiple, subject, channels, time_window, event_types)

%MCP_GET_SUBJECT_EVENTS Returns a matrix that contains Hbo data for target 
%subject in each type, channel, time and type repetition.
%
% First we construct the index matrix that contains specific index of Hbo
% data. Then Use the index matrix to find conresponding Hbo data and add
% them into the output matrix.
%
% Chengyu Deng & Benjamin Zinszer 5 may 2017

%% Index matrix handling (We use MCP struct so we don't need distinguish Homer file version anymore)

% Extract oxy data and marks from the MCP struct
oxy_timeser = mcp_multiple(subject).fNIRS_Data.Hb_data.Oxy(:, channels);
marks_vec = mcp_multiple(subject).fNIRS_Data.Onsets_Matrix;

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
        
        %Abandon the offsets
        temp_marks = temp_marks(1:2:end);
        
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
        %Now abandon the offsets
        temp_marks = temp_marks(1:2:end);
        
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
% The output matrix setup(time x channels x type repetition x types)
event_matrix = nan(length(time_window), length(channels), size(marks_mat, 1), length(event_types));

for type_i = 1 : length(event_types)
    
    matched_conditions = cell2mat(arrayfun(@(x) strcmp(event_types{type_i},x.Name), mcp_multiple(subject).Experiment.Conditions,'UniformOutput',false)); 
    event_marks = [mcp_multiple(subject).Experiment.Conditions(matched_conditions).Mark];
    
    for event_j = 1 : length(marks_mat(:, event_marks))
        if marks_mat(event_j,type_i) + time_window(end) <= length(oxy_timeser),
            event_matrix(:,:,event_j,type_i) = oxy_timeser(marks_mat(event_j,type_i)+time_window,:) - ones(length(time_window),1)*oxy_timeser(marks_mat(event_j,type_i)+time_window(1),:);
        else
            event_matrix(:,:,event_j:end,type_i) = NaN;
        end
    end
end




end