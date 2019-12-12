function summarized_MCPA_struct = summarize_MCPA_Struct(summary_function,MCPA_struct, averaging_dimensions, within_sub, mcp_struct)
%SUMMARIZE_MCPA_STRUCT Convert an MCPA struct of windowed data to
%multivariate patterns using any summarizing function, such as nanmean.
% Fold_type_functions: can be 

if isempty(averaging_dimensions), averaging_dimensions = {'time'}; end

if ischar(summary_function)
    summary_function = str2func(summary_function);
elseif iscell(summary_function)
    summary_function = str2func(summary_function{:});
end

%% get session dimension if this is within subject decoding
if within_sub
    dat_to_summ = add_session_dimension(mcp_struct,MCPA_struct);
else
    dat_to_summ = MCPA_struct.patterns;
end


%% Get summarized:

for dimension = 1:length(averaging_dimensions)
    if dimension == 1
        if strcmp(averaging_dimensions{dimension},'instance')
           try
                summarized_MCPA_struct_pattern = summary_function(dat_to_summ,4);
           catch
                error('Failed to successfully summarize MCPA_struct to a desired form.');
           end      
        elseif strcmp(averaging_dimensions{dimension},'time')
            try
                summarized_MCPA_struct_pattern = summary_function(dat_to_summ,1);
            catch
                error('Failed to successfully summarize MCPA_struct to a desired form.');
            end
        end
    else
        if strcmp(averaging_dimensions{dimension},'instance')
           try
                 summarized_MCPA_struct_pattern = summary_function(summarized_MCPA_struct_pattern,4);
            catch
                error('Failed to successfully summarize MCPA_struct to a desired form.');
            end      
        elseif strcmp(averaging_dimensions{dimension},'time')
            try
                summarized_MCPA_struct_pattern = summary_function(summarized_MCPA_struct_pattern,1);
            catch
                error('Failed to successfully summarize MCPA_struct to a desired form.');
            end
        end
    end
end


summarized_MCPA_struct_pattern = squeeze(summarized_MCPA_struct_pattern);


%% Return the new struct
try
    summarized_MCPA_struct = MCPA_struct;
    summarized_MCPA_struct.created = datestr(now);
    summarized_MCPA_struct.summarizing_function = summary_function;
    summarized_MCPA_struct.patterns = summarized_MCPA_struct_pattern;
    %fprintf(' Done.\n');
    
catch
    error('Failed to create the new struct')
end



end







