function summarized_MCPA_struct = summarize_MCPA_Struct(summary_function,MCPA_struct, summarize_dimensions)
%SUMMARIZE_MCPA_STRUCT Convert an MCPA struct of windowed data to
%multivariate patterns using any summarizing function, such as nanmean.


%% Convert the inputs to correct format
% Convert summary function into a function handle
if ischar(summary_function)
    summary_function = str2func(summary_function);
elseif iscell(summary_function)
    summary_function = str2func(summary_function{:});
elseif isempty(summary_function)
    warning('No summarizing function specified. Reshaping dimensions only.')
end

% Convert summarize_dimensions into cell array
if ischar(summarize_dimensions)
    warning('Please provided cell array of dimensions to summarize. Automatically converting...')
    summarize_dimensions = {summarize_dimensions};
elseif isnumeric(summarize_dimensions)
    warning('Please provided cell array of dimensions to summarize. Looking up dimension names in MCPA''s dimensions field...')
    if ~isfield(MCPA_struct,'dimensions')
        warning('No ''dimension'' field found in MCPA struct, so no dimension labels available. Assigning numeric labels.')
        MCPA_struct.dimensions = cellfun(@(x) strcat('PatternDim',x), cellfun(@num2str, num2cell(1:ndims(MCPA_struct.patterns)),'UniformOutput',false),'UniformOutput',false);
    end
    summarize_dimensions = MCPA_struct.dimensions(summarize_dimensions);
    dim_str = join(summarize_dimensions,', ');
    fprintf('Converted numeric dimension indices to: %s', dim_str{:})
end

%% Create copies of dimensions and of pattern matrix
dimension_labels = MCPA_struct.dimensions;
pattern_matrix = MCPA_struct.patterns;

%% Get summarized
% Step through the dimensions listed from left to right. If a dimension
% contains the concatenation operator (X), then reshape the pattern matrix
% and rename the dimension(s). Otherwise, perform the summarizing function
% on the dimension.

%% new approach
for curr_dim = 1:length(summarize_dimensions)
    
    % Perform concatenations
    concat = strsplit(summarize_dimensions{curr_dim}, 'X');
    if length(concat) == 2
        try
            first_dim_to_concat = find(strcmp(dimension_labels,concat{1}));
            second_dim_to_concat = find(strcmp(dimension_labels,concat{2}));
            [pattern_matrix, new_pattern_matrix_dimensions] = concatenate_dimensions(first_dim_to_concat,second_dim_to_concat,pattern_matrix);
            dimension_labels = new_pattern_matrix_dimensions;
            %summarize_dimensions{curr_dim} = strjoin(concat,'+');
        catch concat_error
            warning(concat_error.message)
            error('Failed to concatenate dimensions: %s', summarize_dimensions{curr_dim})
        end
    elseif ~isempty(summary_function)
        try
            % Perform summary function over the specified dimension
            dim_to_summarize = find(strcmp(dimension_labels, summarize_dimensions{curr_dim}));
            pattern_matrix = summary_function(pattern_matrix,dim_to_summarize);
            
            % Then keep track of remaining dimensions in list
            dimension_labels{dim_to_summarize} = [];
        catch averaging_error
            warning(averaging_error.message);
            error('Failed to successfully summarize dimension: %s', summarize_dimensions{curr_dim});
        end
    end
end
dimension_labels = dimension_labels(~cellfun('isempty',dimension_labels));
summarized_MCPA_struct_pattern = squeeze(pattern_matrix);

%% old approach
% try
%     for dim = 1:length(summarize_dimensions)
%         concat = strsplit(summarize_dimensions{dim}, 'X');
%         if length(concat) == 2
%             first_dim_to_concat = find(strcmp(dimension_labels,concat{1}));
%             second_dim_to_concat = find(strcmp(dimension_labels,concat{2}));
%             [pattern_matrix, new_pattern_matrix_dimensions] = concatenate_dimensions(first_dim_to_concat,second_dim_to_concat,pattern_matrix);
%             dimension_labels = new_pattern_matrix_dimensions;
%         end
%     end
% catch dimension_error
%     warning(dimension_error.message);
%     error('Failed to concatenate these dimensions.');
% end
%
% % Then average over the left over dimensions
% try
%     for dim = 1:length(summarize_dimensions)
%         % skip over the concatenate arguments
%         if length(strsplit(summarize_dimensions{dim}, 'X')) == 2
%             continue;
%         end
%
%         % average over the specified dimension
%         dim_to_average = find(strcmp(dimension_labels, summarize_dimensions{dim}));
%         pattern_matrix = summary_function(pattern_matrix,dim_to_average);
%
%         % Then keep track of remaining dimensions in list
%         dimension_labels{dim_to_average} = [];
%     end
%
%     % clean up dimension cell array
%     dimension_labels = dimension_labels(~cellfun('isempty',dimension_labels));
% catch averaging_error
%     warning(averaging_error.message);
%     error('Failed to successfully summarize MCPA_struct to a desired form. Cannot average over these dimensions');
% end
%
% % Then get rid of all dimensions of length 1
% summarized_MCPA_struct_pattern = squeeze(pattern_matrix);


%% Return the new struct
try
    summarized_MCPA_struct = MCPA_struct;
    summarized_MCPA_struct.created = datestr(now);
    summarized_MCPA_struct.summarizing_function = summary_function;
    summarized_MCPA_struct.patterns = summarized_MCPA_struct_pattern;
    summarized_MCPA_struct.dimensions = dimension_labels;
    fprintf(' Done.\n');
    
catch
    error('Failed to create the new struct');
end



end

