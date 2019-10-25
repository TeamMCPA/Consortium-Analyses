function [mapped_sets,idx_map] = map_values(p, sets)
%% After selecting a subset of channels, map their new indices to what channels they correspond to
% Returns 
% Arguments:
% p: our input struct that results from parse_inputs
% sets: our chosen subsets from find_sets

idx_map = [p.Results.incl_channels; 1:length(p.Results.incl_channels)]';

mapped_sets = zeros(size(sets));
for i = 1:size(sets,1)
    for j = 1:size(sets,2)
        loc = find(idx_map(:,1)==sets(i,j));
        mapped_sets(i,j) = loc;
    end
end
end




