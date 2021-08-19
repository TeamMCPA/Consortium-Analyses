function new_mcp_file = MCP_delete_stimuli(mcp_file,old_label_name,save_flag)
%MCP_DELETE_STIMULI removes stimuli in an MCP struct matching a name or
%trigger number.
% This is useful after the default naming, where integer marks get turned
% into names (e.g., 3 -> 'm3'). The function can be applied to rename
% multiple conditions at once (e.g., 'm2' and 'm3' -> 'cats')
%
% mcp_file may be either a Matlab file (*.mcp extension) or an MCP struct
% in the current workspace.
%
%
% Chengyu Deng & Benjamin Zinszer 7 june 2017
% revised bdz 19 oct 2018

%% Open the old MCP file
if isstruct(mcp_file)
    old_mcp_struct = mcp_file;
else
    [mcpdir, mcpfile, ~] = fileparts(mcp_file);
    old_mcp_struct = load([mcpdir mcpfile '.mcp'],'-mat');
end

%% Recursion for multiple subjects (length MCP > 1)
% If there are multiple subjects, recurse the function for each one without
% saving and then save (if flagged) at the end
if length(old_mcp_struct) > 1
    for subj = 1:length(old_mcp_struct)
        new_mcp_file(subj) = MCP_delete_stimuli(old_mcp_struct(subj),old_label_name,0);
    end
    if save_flag
        save([mcpdir mcpfile '_r.mcp'],'-struct','new_mcp_file')
    end
    return
end

%% Extract the onsets matrix for the condition that will be replaced
% old_cond_index - logical for elements of Conditions that correspond to the old condition.
% old_onset_col - integer array for columns of Onset_Matrix that correspond to the old condition.

% Get the old onset times for any condition(s) that match old name/mark
if isnumeric(old_label_name) % condition where Mark is specified
    % If old_label_name is numeric, check Mark instead of Name
    old_cond_index = arrayfun(@(x)(x.Mark==old_label_name),old_mcp_struct.Experiment.Conditions);
else % otherwise assume Name is specified
    % Otherwise look in the Name field
    old_cond_index = arrayfun(@(x)(strcmp(x.Name,old_label_name)),old_mcp_struct.Experiment.Conditions);
end
old_onset_mark = [old_mcp_struct.Experiment.Conditions(old_cond_index).Mark];


%% Mask out the column in onsets and remove Condition entry
% This will delete the Condition and leave the column behind, but that is
% expedient for now because otherwise we'd have to update every Mark to
% match the new column numbers in the onset matrix
new_mcp_file = old_mcp_struct;
new_mcp_file.Experiment.Conditions(old_cond_index)=[];
new_mcp_file.fNIRS_Data.Onsets_Matrix(:,old_mcp_struct.Experiment.Conditions(old_cond_index).Mark)=0;

%% If the save_flag is true, write the data out.
if save_flag, save([mcpdir mcpfile '_r.mcp'],'-struct','new_mcp_file'); end
