
function new_patterns = concatenate_dimensions(patterns, concat_dims, varargin)
% reshape a matrix (patterns) depending on concat_dims
% concat_dims should be a numeric value indicating which dimensions to concatenate over
% ex. if we had the had a matrix of 4 dimensiosn (1 2 3 4) [2 3] will concatenate dimensions 2 and 3 such that our 
% matrix is now three dimensions (1 [2 3] 4) 


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

% Keep dimensions that are dimensions preserved from the original pattern
% matrix and not concatenated.
keep_dims = 1:ndims(patterns);
% Remove the dimensions that will be concatenated from keep_dims list
[~, drop_dims, ~] = intersect(keep_dims,concat_dims);
keep_dims(drop_dims) = [];

% Permute the patterns matrix to put the concatenated dimensions at the end
patterns = permute(patterns,[keep_dims, concat_dims]);

% Reshape the patterns matrix to stack up the concatenated dimensions
new_patterns = reshape(patterns,[old_pattern_size(keep_dims),new_dimension_size]);

% put back in proper order
rearranged_dims = [keep_dims, concat_dims]; % this is what the dimensions currently look like - all concatenated dims are at the end

amount_removed_dims = length(old_pattern_size) - ndims(new_patterns); % find how many dimensions we removed by looking at difference between the dimensions

rearranged_dims((end-amount_removed_dims+1):end) = []; % rearranged_dims still has whatever dimensions we removed, so we remove those from our record keeping 

if max(rearranged_dims) > length(rearranged_dims)
    [~, idx] = max(rearranged_dims); % the maximum dimension is no longer the true max. meaning if we had 6 dimensions and went down to 5, we're still counting that 6th dimension as 6 and not 5
    rearranged_dims(idx) = rearranged_dims(idx) - amount_removed_dims; % so we subtract
end

mapping = [rearranged_dims; 1:length(rearranged_dims)]; % this creates a map between what we think of as the rearranged dimensions and how matlab sees the rearranged dimensions (i.e. we see it as 1 2 3 5 4, matlab sees 1 2 3 4 5)
[temp, order] = sort(mapping(1,:)); % sort out the rearranged dims so they are in increasing order
sorted_dims = mapping(:,order); % use this rearrangement to tell matlab how to reshape our dimensions     

new_patterns = permute(new_patterns, sorted_dims(2,:)); % reshape
    
end
