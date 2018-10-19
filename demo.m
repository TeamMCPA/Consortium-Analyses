% Get a list of NIRS files
nirs_files = arrayfun(@(x) x.name,dir('*.nirs'),'UniformOutput',false);

% Filter out any preprocessed files (*_proc.nirs)
nirs_files_raw = nirs_files(cellfun(@(x) isempty(x),strfind(nirs_files,'_proc.nirs')));
nirs_files_proc = setdiff(nirs_files,nirs_files_raw);

% Preprocess the raw nirs files and get a list of preprocessed filenames
if length(nirs_files_proc) < length(nirs_files_raw),
    nirs_files_proc_today = preproc_hmr(nirs_files_raw,'Pton_standard_20180924.mat');
    nirs_files_proc = [nirs_files_proc; nirs_files_proc_today];
end

% Build an MCP struct to describe all of the preprocessed data
MCP_data = build_MCP(...
    {nirs_files_proc(1:4),nirs_files_proc(5:8)},...     % the filenames
    {'subject02','subject04'},...                       % subject IDs
    {'shimadzu139','shimadzu139'},...                   % probe IDs
    's');                                               % field for stims


%cond_names = {'baby', 'book', 'bottle', 'cat', 'dog', 'hand', 'shoe', 'spoon'};
cond_names = {'anim', 'inan', 'inan', 'anim', 'anim', 'anim', 'inan', 'inan'};

for old_mark = 1:8
    MCP_data = MCP_relabel_stimuli(MCP_data,old_mark,cond_names{old_mark},0);
end

% Save the compiled dataset
save(['MCP_data_' date '.mcp'],'MCP_data');

% Extract events into MCPA struct
MCPA_data = MCP_to_MCPA(MCP_data,[1:2],[1:139],[0:60]);
save(['MCPA_data_' date '.mat'],'MCPA_data');
