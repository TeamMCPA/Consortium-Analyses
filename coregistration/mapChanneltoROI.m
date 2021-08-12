function transformation_mat = mapChanneltoROI(n_chan, brain_atlas, channel_locs, areas, use_proportional)
%% co-register channels to a region of interest
% n_chan = number of channels on nirs cap, not number included in analysis
% n_areas = number of regions of interest
% brain_atlas = file path to find the TDdatabase
% channel_locs = channel_locs file that contains where the probes sat for a single
% participant's session
% use_proportional = true/false on whether to use proportional channel assignment

%% initialize variables
n_areas = length(areas);
transformation_mat = zeros(n_chan, n_areas);
brain_map = brain_atlas;

%% then find distance for every channel
channel_locations = channel_locs.R.ch.xyzC';


if use_proportional
    rad = 15; 
    for chan = 1:n_chan    
        for area = 1:length(areas)
            eval(['areaMNI = wholeMaskMNIAll.',areas{area},';']);
            transformation_mat(chan, area) = sum(sqrt(sum((areaMNI - channel_locations(chan,:)).^2,2)) <= rad);
        end
        transformation_mat(chan, :) = transformation_mat(chan, :)./sum(transformation_mat(chan, :));
    end
else
    for chan = 1:n_chan
    areaDist = nan(1,n_areas);
    for area = 1:n_areas
        eval(['areaMNI = brain_map.',areas{area},';']);
        if isempty(areaMNI)
            continue;
        end
        areaDist(area) = min(sqrt(sum((areaMNI - repmat(channel_locations(chan,:),size(areaMNI,1),1)).^2,2)));
    end
    %areaDist = areaDist(areaDist ~= 0);
    [val, chan_area] = min(areaDist);
    transformation_mat(chan, chan_area) = 1;
end

end
