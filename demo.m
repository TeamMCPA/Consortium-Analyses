% Get a list of NIRS files
nirs_files = arrayfun(@(x) x.name,dir('*.nirs'),'UniformOutput',false);

% Filter out any preprocessed files (*_proc.nirs)
nirs_files_raw = nirs_files(cellfun(@(x) isempty(x),strfind(nirs_files,'_proc.nirs')));
nirs_files_proc = setdiff(nirs_files,nirs_files_raw);

% Generate a list of files we can use in the MCP builder
[nirs_files_forMCP, unique_subjs] = prep_nirsfiles_mcp(nirs_files_proc,'_','subject');

% Build an MCP struct to describe all of the preprocessed data
MCP_data = build_MCP(...
    nirs_files_forMCP,...                               % the filenames
    unique_subjs,...                                       % subject IDs
    repmat({'shimadzu139'},length(unique_subjs),1),...   % probe IDs
    's');                                               % field for stims

%cond_names = {'baby', 'book', 'bottle', 'cat', 'dog', 'hand', 'shoe', 'spoon'};
cond_names = {'anim', 'inan', 'inan', 'anim', 'anim', 'anim', 'inan', 'inan'};

for old_mark = 1:8
    MCP_data = MCP_relabel_stimuli(MCP_data,old_mark,cond_names{old_mark},0);
end

% Save the compiled dataset
MCP_file_name = ['MCP_data_' date '.mcp'];
save(MCP_file_name,'MCP_data');

% Perform individual event classification (MCPA method)
eventResults = nfold_classify_IndividualEvents(MCP_data, ...
    'cond1','anim','cond2','inan', ...
    'time_window',[2,6], ...
    'summary_handle',@nanmean,'test_handle',@mcpa_classify, ...
    'setsize',139);

% Extract events into MCPA struct and summarize by taking average over time
MCPA_data = MCP_to_MCPA(MCP_data,[],[],[2:6]);
save(['MCPA_data_' date '.mat'],'MCPA_data');
MCPA_data_summarized = summarize_MCPA_Struct(@nanmean,MCPA_data,1);

% Between Subjects decoding
Between_Subj_Result = leave_one_Ss_out_classifyAverages(MCPA_data_summarized,1,2);
