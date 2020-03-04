function transformation_mat = mapChanneltoROI(n_chan, n_features, subj_idx, session_idx)
%% create list of area names
load('TDdatabase.mat');
areas = {};
for area = 1:47
    area_name = ['brodmann_area_', int2str(area)];
    areas{end+1} = area_name;
end

%% load the POS file
pathname = ['/Users/annaph/Desktop/decoding_POS_data/dec',num2str(subj_idx,'%02.0f'),'/decoding_',num2str(session_idx),'/'];
load([pathname 'POS.mat']);

%% initialize the transformation matrix
transformation_mat = zeros(n_chan, n_features);

%% then find distance for every channel
channel_locations = R.ch.xyzC';

for chan = 1:n_chan
    for area = 1:length(areas)
        eval(['areaMNI = wholeMaskMNIAll.',areas{area},';']);
        if isempty(areaMNI)
            continue;
        end
        areaDist(area) = min(sqrt(sum((areaMNI - repmat(channel_locations(chan,:),size(areaMNI,1),1)).^2,2)));
    end
    areaDist = areaDist(areaDist ~= 0);
    [val, chan_area] = min(areaDist);
    transformation_mat(chan, chan_area) = 1;
end

%% then weight by how many channels are in each area
weight = sum(transformation_mat,1);
transformation_mat = transformation_mat ./ weight;
transformation_mat(isnan(transformation_mat))=0;

end

