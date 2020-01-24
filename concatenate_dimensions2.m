function [new_patterns, new_dimensions] = concatenate_dimensions2(patterns, concat_dims, varargin)

% In case somebody enters the dims as a list of arguments instead of a
% single vector
if nargin>2
    concat_dims = [concat_dims, cell2mat(varargin)];
end

% Determine the current dimensions (used in reshape later)
old_pattern_size = size(patterns);

% Size of the new dimension will be the product of the sizes of all
% of the concatenated dimensions
new_dimension_size = prod(old_pattern_size(concat_dims));

% Keep dimensions are the dimensions preserved from the original pattern
% matrix and not concatenated.
keep_dims = 1:ndims(patterns);
% Remove the dimensions that will be concatenated
[~, drop_dims, ~] = intersect(keep_dims,concat_dims);
keep_dims(drop_dims) = [];

% Permute the patterns matrix to put the concatenated dimensions at the end
patterns = permute(patterns,[keep_dims, concat_dims]);

% Reshape the patterns matrix to stack up the concatenated dimensions
new_patterns = reshape(patterns,[old_pattern_size(keep_dims),new_dimension_size]);
new_dimensions = [keep_dims, concat_dims];