function summarized_MCPA_struct = summarize_MCPA_Struct(function_name,MCPA_struct, dimension)
%SUMMARIZE_MCPA_STRUCT Convert an MCPA struct of windowed data to
%multivariate patterns using any summarizing function, such as nanmean.
%
%The function is called with the following arguments:
%function name: either a function handle (e.g., @nanmean) or a character
%array containing the function name (e.g., 'nanmean').
%MCPA_struct: An MCPA struct containing a pattern matrix.
%dimension: Which dimension of the pattern matrix to operate over.
%Available dimensions are (time x event_type x channel x subject). If the
%dimension field is left blank, it will default to the first dimension.
%
%The function will return a new MCPA struct whose field patterns are
%summarized by the specified method.
%
% Chengyu Deng & Benjamin Zinszer 5 may 2017
% revised bdz 26 oct 2018

if ~exist('dimension','var'), dimension = 1; end

if ischar(function_name)
    function_name = str2func(function_name);
elseif iscell(function_name)
    function_name = str2func(function_name{:});
end

%% Get summarized:
try
    summarized_MCPA_struct_pattern = squeeze(function_name(MCPA_struct.patterns,dimension));
catch
    error('Failed to successfully summarize MCPA_struct to a desired form.');
end

%% Return the new struct
try
    summarized_MCPA_struct = MCPA_struct;
    summarized_MCPA_struct.created = datestr(now);
    summarized_MCPA_struct.summarizing_function = function_name;
    summarized_MCPA_struct.patterns = summarized_MCPA_struct_pattern;
    fprintf(' Done.\n');
    
catch
    error('Failed to create the new struct')
end



end







