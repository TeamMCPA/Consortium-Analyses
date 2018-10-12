% Get a list of NIRS files
nirs_files = arrayfun(@(x) x.name,dir('*.nirs'),'UniformOutput',false);

% Filter out any preprocessed files (*_proc.nirs)
nirs_files_raw = nirs_files(cellfun(@(x) isempty(x),strfind(nirs_files,'_proc.nirs')));

% Preprocess the nirs files and get a list of preprocessed filenames
nirs_files_proc = preproc_hmr(nirs_files_raw,'Pton_standard_20180924.mat');

% Build an MCP struct to describe all of the preprocessed data
MCP_data = build_MCP(...
    {nirs_files_proc(1:4),nirs_files_proc(5:8)},...     % the filenames
    {'subject02','subject04'},...                       % subject IDs
    {'shimadzu139','shimadzu139'},...                   % probe IDs
    's');                                               % field for stims


