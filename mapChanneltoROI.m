function transformation_mat = mapChanneltoROI(n_chan, database, pos, areas)
%% create a matrix of 0's and 1's indicating which anatomical region a channel corresponds to
% n_chan = number of channels on nirs cap, not number included in analysis
% database = atlas containing MNI coordinates for Brodmann's areas
% pos = POS file that contains where the probes sat for a single
% areas = the fieldnames in the database struct that correspond to the Brodmann's areas we want to co-register to


%% initialize the transformation matrix
n_areas = length(areas);
transformation_mat = zeros(n_chan, n_areas);

%% initialize other variables
brain_map = database;

%% then find distance for every channel
channel_locations = pos.R.ch.xyzC';

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







