function [nirs_files_sorted, subject_ids] = prep_nirsfiles_mcp(nirs_files_list,delim_string,subj_string)

% If nirs_files_list is a string, assume it's a searchstring for dir()
if isstring(nirs_files_list) || ischar(nirs_files_list)
    nirs_files_list = dir(nirs_files_list);
end

% If nirs_files_list is a struct (like the output of dir search) then convert
if isstruct(nirs_files_list)
    nirs_files_list = arrayfun(@(x) [x.folder filesep x.name],nirs_files_list,'UniformOutput',false);
end

% Determine which subjects each file belongs with
[~, file_names, ~] = cellfun(@fileparts, nirs_files_list,'UniformOutput',false);
file_names_parsed = cellfun(@(x) strsplit(x,delim_string),file_names,'UniformOutput',false);
subj_list = cellfun(@(x) x{contains(x,subj_string)},file_names_parsed,'UniformOutput',false);


% Create file array
subject_ids = unique(subj_list);
nirs_files_sorted = cell(length(subject_ids),1);
for s_idx = 1:length(subject_ids)
    nirs_files_sorted{s_idx} = nirs_files_list(strcmp(subj_list,subject_ids(s_idx)));
end
    
