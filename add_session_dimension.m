function session_mcpa_pattern = add_session_dimension(mcp_data,mcpa_struct)
% start dimensions: time x condition x channel x instance x subject

%% first get max amount of trials in a session
    % get max sessions while you do this

numSessions = nan(length(mcp_data),1); % initialize empty array to store how many sessions
for subject = 1:length(mcp_data)
    numSessions(subject) = length(mcp_data(subject).Experiment.Runs);
end
max_sessions = max(numSessions);


numTrials = nan(length(mcp_data),max(numSessions)); % initialize empty matrix to store max trials in each session
for subject = 1:length(mcp_data)
    for session = 1:length(mcp_data(subject).Experiment.Runs)
        session_idx = mcp_data(subject).Experiment.Runs(session).Index'; % get the indices for that sessions onsets matrix
        stims=(mcp_data(subject).fNIRS_Data.Onsets_Matrix(session_idx,:) * [1;2;3;4;5;6;7;8]); % get the stims in form of their mark
        locs = find(stims);
        numTrials(subject, session) = length(locs);   
    end
end
max_trials_in_session = max(max(numTrials));        


%% then need to know what session each item is in 
session_data = nan(max_trials_in_session, max_sessions, length(mcp_data)); % empty matrix to store data in 

num_marks = length(mcp_data(1).Experiment.Conditions);
for subject = 1:length(mcp_data) % for each subject
    for session = 1:length(mcp_data(subject).Experiment.Runs) % for each session
        session_idx = mcp_data(subject).Experiment.Runs(session).Index'; % get the indices for that sessions onsets matrix
        stims=(mcp_data(subject).fNIRS_Data.Onsets_Matrix(session_idx,:) * [1:num_marks]'); % get the stims in form of their mark
        locs = find(stims); % get stim locations
        nans_to_pad = nan(length(session_data) - length(stims(locs)),1); % need to pad in nans since not all sessions have same amount of trials
        that_sessions_data = vertcat(stims(locs), nans_to_pad); 
        session_data(:,session,subject) = that_sessions_data;
    end
end

%% now need counts for each event in each session
% this will let us know what indicies to 'skim' events off from in the 

event_count_in_session = zeros(num_marks,max_sessions,length(mcp_data)); % create empty matrix of ????
for sub = 1:length(mcp_data) % for each subject
    for session = 1:length(mcp_data(sub).Experiment.Runs) % in each of their sessions
        tbl = tabulate(session_data(:,session,sub)); % count the data in that session
        if isempty(tbl) % if that session had no data (since some have fewer sessions)
            event_count_in_session(:,session,sub) = nan(num_marks,1); % fill in 
        else % if that session does have data
            if length(tbl) < num_marks % see if that data has less than 8 categories so that we can pad the vector with nans
                nans_to_pad = zeros(num_marks-size(tbl,1), 1);
                eventCount = vertcat(tbl(:,2),nans_to_pad);
                
                event_count_in_session(:,session,sub) = eventCount;
            else % otherwise we can just do straight forward assignment
                event_count_in_session(:,session, sub) = tbl(:,2);
            end
        end
    end
end

%% now we assign values to each fold

dim1 = size(mcpa_struct.patterns,1); % time samples
dim2 = size(mcpa_struct.patterns,2); % conditions
dim3 = size(mcpa_struct.patterns,3); % channels
dim4 = size(mcpa_struct.patterns,4); % instances
dim5 = max_sessions;                 % sessions
dim6 = size(mcpa_struct.patterns,5); % participant

session_mcpa_pattern = nan(dim1,dim2,dim3,dim4,dim5,dim6);


for subject = 1:length(mcp_data)
    for session = 1:max_sessions
        for cond = 1:num_marks
            if ~all(isnan(event_count_in_session(:,session,sub)))
                instance_idx = event_count_in_session(cond,1:session,sub);
                max_idx = sum(instance_idx);
                min_idx = max_idx - instance_idx(end)+1;
                
                pattern_mat = mcpa_struct.patterns(1:dim1,cond,1:dim3,min_idx:max_idx,sub);
                temp_mat = pattern_mat(:,:,:,:,1,1);
                session_mcpa_pattern(:,cond,:, min_idx:max_idx, session,subject) = temp_mat;
                % need to figure out how to do assignment here
                
            end
        end
    end
end

% % now need to begin selecting values for each fold
% session_struct = struct;
% for sub = 1:length(mcp_data)
%     for session = 1:4
%         for cond = 1:8
%             if ~all(isnan(y(:,session,sub)))
% 
%                 instance_idx = y(cond,1:session,sub);
%                 max_idx = sum(instance_idx);
%                 min_idx = max_idx - instance_idx(end)+1;
% 
%                 session_struct(sub).subj(session).session(cond).condition.pattern_dat = mcpa_struct.patterns(:,cond,:,min_idx:max_idx,sub);
%             end
%         end
%     end
% end
    

% final dimensions: time x condition x channel x repetition x session x subject


end
