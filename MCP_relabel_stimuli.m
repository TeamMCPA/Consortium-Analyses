function new_mcp_file = MCP_relabel_stimuli(mcp_file,old_label_name,new_labels,save_flag)
%MCP_RELABEL_STIMULI change names of stimuli in an MCP struct
% The function is called with the following arguments:
% MCP_relabel_stimuli(mcp_file,old_label_name,new_labels,save_flag)
%
% This is useful after the default naming, where integer marks get turned
% into names (e.g., 3 -> 'm3'). The function can be applied to rename
% multiple conditions at once (e.g., 'm2' and 'm3' -> 'cats')
%
% mcp_file: may be either a Matlab file (*.mcp extension) or an MCP struct
% in the current workspace.
% 
% old_label_name: may be either numeric trigger number or a char-type
% containing the previous name (e.g., 'm3' or 'cats')
%
% new_labels: may be either a char-type containing the new name (e.g.,
% 'felines') or a cell array of replacement names (e.g.,
% {'tabby','calico','tiger','garfield'}). However the cell array must have
% the same length as the number of instances of the old_label_name found in
% the mcp_file, or it will throw an error.
%
% save_flag: logical, set whether to write out a new *.mcp file. The file
% with have the same name as the mcp_file input, but with _r suffix.
%
%
% Chengyu Deng & Benjamin Zinszer 7 june 2017
% revised bdz 19 oct 2018

%% Clean up / process inputs
% If necessary, convert the new_labels variable into a cell array
if isnumeric(new_labels)
    new_labels = cellstr(num2str(new_labels));
elseif ischar(new_labels)
    new_labels = cellstr(new_labels);
end

% Open the old MCP file
if isstruct(mcp_file)
    old_mcp_struct = mcp_file;
else
    [mcpdir, mcpfile, ~] = fileparts(mcp_file);
    %fprintf('Loading file: %s\n\n',mcpfile);
    old_mcp_struct = load([mcpdir mcpfile '.mcp'],'-mat');
end

%% Recursion for multiple subjects (length MCP > 1)
% If there are multiple subjects, recurse the function for each one without
% saving and then save (if flagged) at the end

% Quickly create char version of old label for convenient screen display
if isnumeric(old_label_name)
    old_label_char = num2str(old_label_name);
else
    old_label_char = old_label_name; 
end

if length(old_mcp_struct) > 1
    for subj = 1:length(old_mcp_struct)
        disp(['Updating labels "' old_label_char '" for: ' old_mcp_struct(subj).Subject.Subject_ID] );
        new_mcp_file(subj) = MCP_relabel_stimuli(old_mcp_struct(subj),old_label_name,new_labels,0);
    end
    if save_flag
        disp(['Saving new file under name: ' mcpfile '_r.mcp']);
        disp([]);
        save([mcpdir mcpfile '_r.mcp'],'-struct','new_mcp_file')
    end
    disp([]);
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
old_onset_col = [old_mcp_struct.Experiment.Conditions(old_cond_index).Mark];
old_cond_onsets = old_mcp_struct.fNIRS_Data.Onsets_Matrix(:,old_onset_col);

%% Compare the new labels and the old label
% If the number of new labels is not equal to the number of onsets in the
% condition they are replacing, two possible cases may occur: (1) all of the
% old label items are being replaced with a single new label (e.g., all m1
% become "dog") or (2) the trials are supposed to be relabeled one-by-one,
% but there was a mistake in the number of replacements specified. In the
% first case, we can just do a big batch-replace. In the second case, it's
% just an error and needs to be bounced back.

if length(new_labels) ~= sum(old_cond_onsets>0)
    % Case (1): Prepare to replace all instances of old with the new label.
    if length(new_labels) == 1
        disp(['Only one new label provided for ' num2str(sum(old_cond_onsets)) ' onsets. Applying label to all onsets.']);
        new_mcp_file = old_mcp_struct;
        new_mcp_file.Experiment.Conditions(old_cond_index).Name = new_labels{:};
        return
        % This was the old way of doing it, less efficient but reuses the
        % code below for individual trial labeling. It also creates new
        % columns even though the new is exact duplicate of the old.
        %new_labels = repmat(new_labels,max(sum(old_cond_onsets)),1);
    elseif length(new_labels)==1 && sum(old_cond_onsets)==1
        % Case (2): In the occasional instance that there is only one onset, it
        % won't meet the criteria above, but we still just want to do a label
        % swap in the Conditions rather than adding a new column.
        disp('Only one onset found for this label. Replacing');
        new_mcp_file = old_mcp_struct;
        new_mcp_file.Experiment.Conditions(old_cond_index).Name = new_labels{:};
        return
       
    % Case (3): Throw an error and quit.
    else
        new_mcp_file = old_mcp_struct;
        error('Provided %g new labels to replace %g items in existing condition! These values must be equal.',length(new_labels),sum(old_cond_onsets));
        return
    end
    
    
    
end

% Turn the list of labels into a set of integers which can be acted upon.
num_existing_conds = length(old_mcp_struct.Experiment.Conditions);
[ unique_new, unique_integers, new_integer_labels ] = unique(new_labels);
num_new_conds = length(unique_integers);

% Fill out a new set of onsets, first as integers in a single vector.
new_cond_onsets = old_cond_onsets;
new_cond_onsets(new_cond_onsets==1) = new_integer_labels;
% Second, as a matrix with a column for each integer-label and logicals for
% the onsets of that integer-label in the time series
new_cond_onsets = repmat(new_cond_onsets,1,num_new_conds) == repmat(new_integer_labels(unique_integers)',size(new_cond_onsets));

% Copy the old MCP file
new_mcp_file = old_mcp_struct;

% Append the new onsets matrix to the right edge of the old onsets matrix
new_mcp_file.fNIRS_Data.Onsets_Matrix = [old_mcp_struct.fNIRS_Data.Onsets_Matrix new_cond_onsets];

for new_cond = 1:length(unique_new)
    % Append the new condition names to the list of conditions
    new_mcp_file.Experiment.Conditions(num_existing_conds+new_cond).Name = unique_new{new_cond};
    % Assign a new mark number that refers to the column of Onsets_Matrix
    new_mcp_file.Experiment.Conditions(num_existing_conds+new_cond).Mark = size(old_mcp_struct.fNIRS_Data.Onsets_Matrix,2)+new_cond;
end

% After conversion is finished, delete the old condition. It turns out this
% is kind of necessary because otherwise it gums up the RSA matrix later
%new_mcp_file.Experiment.Conditions(old_cond_index)=[];
%new_mcp_file.fNIRS_Data.Onsets_Matrix(:,old_mcp_struct.Experiment.Conditions(old_cond_index).Mark)=[];
% Turns out this was a bad idea because is scrambles up the correspondance
% between Mark # and the column of Onsets_Matrix

%% If the save_flag is true, write the data out.
if save_flag, save([mcpdir mcpfile '_r.mcp'],'-struct','new_mcp_file'); end