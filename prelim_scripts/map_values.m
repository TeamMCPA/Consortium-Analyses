function [mapped_sets,idx_map] = map_values(p, sets)

idx_map = [p.Results.incl_channels; 1:length(p.Results.incl_channels)]';

mapped_sets = zeros(size(sets));
for i = 1:size(sets,1)
    for j = 1:size(sets,2)
        loc = find(idx_map(:,1)==sets(i,j));
        mapped_sets(i,j) = loc;
    end
end
end




