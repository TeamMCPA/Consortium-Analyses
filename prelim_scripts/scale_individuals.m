function mcp_struct = scale_individuals(mcp_struct, input_struct)
%% Scale individual participants' data
% Arguments: 
% MCP data struct containing all pariticpant's data
% input_struct: input that we entered into the folding function with our
% classification parameters

if input_struct.norm_withinSessions % if we want to norm within each session
    for participant = 1:length(mcp_struct) % for each participant
        scaled_data = [];
        for session = 1:length(mcp_struct.Experiment.Runs) % scale each session
            current_dataset_idx = mcp_struct.Experiment.Runs(session).Index'; 
            current_dataset= mcp_struct.fNIRS_Data.Hb_data.Oxy(current_dataset_idx,:);
            session_scaled = input_struct.norm_function(current_dataset,[],[], input_struct);
            scaled_data = vertcat(scaled_data, session_scaled);
        end
        mcp_struct(participant).fNIRS_Data.Hb_data.Oxy = scaled_data; % then recombine
    end
else % if we want to norm across sessions
    for participant = 1:length(mcp_struct) % for each participant, scale their whole Oxy dataset
        scaled_data = input_struct.norm_function(mcp_struct.fNIRS_Data.Hb_data.Oxy,[],[], input_struct);
        mcp_struct(participant).fNIRS_Data.Hb_data.Oxy = scaled_data;
    end
end

end