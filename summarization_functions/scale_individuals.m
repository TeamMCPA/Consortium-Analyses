function mcp_struct = scale_individuals(mcp_struct, input_struct)
%% Scale individual participants' data
% Arguments: 
% MCP data struct containing all pariticpant's data
% input_struct: input that we entered into the folding function with our
% classification parameters

if ~isfield(input_struct, 'scale_withinSessions')
    input_struct.scale_withinSessions = true;
end

%% figure out which data we are scaling
hemo_types = strsplit(input_struct.hemoglobin, '+');
for hemo_type = 1:length(hemo_types)
    current_hemo_type = hemo_types{hemo_type};
    for participant = 1:length(mcp_struct)
        % create empty matrix to store results
        scaled_hemo_data = nan(size(mcp_struct(participant).fNIRS_Data.Hb_data.(current_hemo_type))); 
        for session = 1:length(mcp_struct(participant).Experiment.Runs) % scale each session individually
            % get indices for this dataset
            current_dataset_idx = mcp_struct(participant).Experiment.Runs(session).Index';            
            % scale data
            session_scaled = input_struct.scale_function(mcp_struct(participant).fNIRS_Data.Hb_data.(current_hemo_type)(current_dataset_idx,:), input_struct);            
            % save to new matrix
            scaled_hemo_data(current_dataset_idx,:) = session_scaled;  
        end
        mcp_struct(participant).fNIRS_Data.Hb_data.(current_hemo_type) = scaled_hemo_data;
    end
end
    
end
