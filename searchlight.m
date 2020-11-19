function [chan_list] = searchlight(MCP_data,seed,N)
% Searchlight function for selecting the N closest channels to seed
% channel(s)
% inputs:
% MCP_data: multivariate pattern data  (data, with correct src and detector
% coordinates)
% seed: vector channel number(s) corresponding to the desired center points 
% from which to compute the searchlight. If empty > all channels
% N: number of surrounding channels
% outputs: 
% chan_list: length(seed) x N  matrix storing the N closest channels 
% surrounding each of the length(seed) seed channels
% Alexis Black, Claire Kabdebon, Alice Wang -- November 2020

%if empty then assign it all channels
if isempty(seed)
    seed = MCP_data(1).Experiment.Probe_arrays.Channels;
end

%first compute the channel coordinates : here we postulate that couples(1,:)
%corresponds to channel #1: that should be verified !
%retrieve the Src Det couples 
channels = MCP_data(1).Experiment.Probe_arrays.Channels;
couples = MCP_data(1).Experiment.Probe_arrays.Geometry.MeasList(channels,1:2);

src = couples(:,1);
det = couples(:,2);

src_coord = MCP_data(1).Experiment.Probe_arrays.Geometry.SrcPos(src,:);
det_coord = MCP_data(1).Experiment.Probe_arrays.Geometry.DetPos(det,:);

channels = (src_coord+det_coord)/2;

%initialize chan_list
chan_list = nan(length(seed),N);
%go through seeds
for s = 1:length(seed)
    %retrieve the seed coordinates
    seed_coord = channels(seed(s),:);

    %for each channel compute the distance from the seed
    dist=zeros(size(channels,1),1);
    for ch=1:size(channels,1)
        dist(ch) = norm(channels(ch,:)-seed_coord);
    end

    % sort the distances from smallest to largest
    [~,indices] = sort(dist,'ascend');

    %retrieve the N closest channels to the seed
    chan_list(s,:) = indices(1:N);
end

end