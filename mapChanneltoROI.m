function transformation_mat = mapChanneltoROI(n_chan, database, pos, areas)
% n_chan = number of channels on nirs cap, not number included in analysis
% n_areas = number of Brodmann's areas - eventually we could update this so
% that we can also use other regions of interest, but for now this parameter
% should always be 47
% td_path = file path to find the TDdatabase
% pos = POS file that contains where the probes sat for a single
% participant's session


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







