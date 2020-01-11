function summarized_MCPA_struct = summarize_MCPA_Struct(summary_function,MCPA_struct, summarize_dimensions)
%SUMMARIZE_MCPA_STRUCT Convert an MCPA struct of windowed data to
%multivariate patterns using any summarizing function, such as nanmean.

if ischar(summary_function)
    summary_function = str2func(summary_function);
elseif iscell(summary_function)
    summary_function = str2func(summary_function{:});
end

%% create copy of dimensions and of pattern matrix
dimensions = MCPA_struct.dimensions;
pattern_matrix = MCPA_struct.patterns;

%% Get summarized:

% search for what dimensions to concatenate and keep record of
% new dimensions
try
    for dim = 1:length(summarize_dimensions)
        concat = strsplit(summarize_dimensions{dim}, 'X');
        if length(concat) == 2
            first_dim_to_concat = find(strcmp(dimensions,concat{1}));
            second_dim_to_concat = find(strcmp(dimensions,concat{2}));
            [pattern_matrix, new_pattern_matrix_dimensions] = concatenate_dimensions(first_dim_to_concat,second_dim_to_concat,pattern_matrix);
            dimensions = new_pattern_matrix_dimensions;
        end  
    end
catch dimension_error
    warning(dimension_error.message);
    error('Failed to concatenate these dimensions.');
end

% Then average over the left over dimensions
try
    for dim = 1:length(summarize_dimensions)
        % skip over the concatenate arguments
        if length(strsplit(summarize_dimensions{dim}, 'X')) == 2
            continue;
        end
        
        % average over the specified dimension
        dim_to_average = find(strcmp(dimensions, summarize_dimensions{dim}));
        pattern_matrix = summary_function(pattern_matrix,dim_to_average);
        
        % Then keep track of remaining dimensions in list
        dimensions{dim_to_average} = [];
    end
    
    % clean up dimension cell array
    dimensions = dimensions(~cellfun('isempty',dimensions));
catch averaging_error
    warning(averaging_error.message);
    error('Failed to successfully summarize MCPA_struct to a desired form. Cannot average over these dimensions');
end

% Then get rid of all dimensions of length 1
summarized_MCPA_struct_pattern = squeeze(pattern_matrix);


%% Return the new struct
try
    summarized_MCPA_struct = MCPA_struct;
    summarized_MCPA_struct.created = datestr(now);
    summarized_MCPA_struct.summarizing_function = summary_function;
    summarized_MCPA_struct.patterns = summarized_MCPA_struct_pattern;
    summarized_MCPA_struct.dimensions = dimensions;
    fprintf(' Done.\n');

catch 
    error('Failed to create the new struct');
end



end

