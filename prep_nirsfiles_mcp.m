function [nirs_files_sorted, subject_ids] = prep_nirsfiles_mcp(nirs_files_list,delim_string,subj_string)

% Determine which subjects each file belongs with
file_names_parsed = cellfun(@(x) strsplit(x,delim_string),nirs_files_list,'UniformOutput',false);
subj_list = cellfun(@(x) x{contains(x,subj_string)},file_names_parsed,'UniformOutput',false);


% Create file array
subject_ids = unique(subj_list);
nirs_files_sorted = cell(length(subject_ids),1);
for s_idx = 1:length(subject_ids)
    nirs_files_sorted{s_idx} = nirs_files_list(strcmp(subj_list,subject_ids(s_idx)));
end
    