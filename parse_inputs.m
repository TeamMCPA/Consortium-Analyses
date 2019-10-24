function p = parse_inputs(MCP_struct, varargin)
%% input parser 
% takes input from function that calls it and creates a struct to store
% clasification parameters
% if not input is provided for a particular variable, it creates default values

% Arguments:
% MCP_struct: MCP data structure
% varargin: input to main function

p = inputParser;
addParameter(p,'incl_channels',[1:max(arrayfun(@(x) size(x.fNIRS_Data.Hb_data.Oxy,2),MCP_struct))],@isnumeric);
addParameter(p,'incl_subjects',[1:length(MCP_struct)],@isnumeric);
addParameter(p,'time_window',[2,6],@isnumeric);
addParameter(p,'baseline_window',[-5 0],@isnumeric);
addParameter(p,'conditions',unique(cellstr(char(cellfun(@(x) char(x{:}), arrayfun(@(x) unique({x.Experiment.Conditions.Name},'stable'),MCP_struct, 'UniformOutput',false),'UniformOutput',false))),'stable'),@iscell);
addParameter(p,'summary_handle',@nanmean);
addParameter(p,'setsize',max(arrayfun(@(x) size(x.fNIRS_Data.Hb_data.Oxy,2),MCP_struct)),@isnumeric);
addParameter(p,'max_sets',1000000,@isnumeric);
addParameter(p,'test_handle',@mcpa_classify);
addParameter(p,'opts_struct',[],@isstruct);
addParameter(p,'verbose',true,@islogical);
addParameter(p,'norm_data', false, @islogical);
addParameter(p, 'norm_function', @minMax_scale_0to1);

parse(p,varargin{:});

end
