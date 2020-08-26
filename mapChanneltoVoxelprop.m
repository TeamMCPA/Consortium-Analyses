function transformation_mat = mapChanneltoVoxelprop(n_chan, n_areas, nii_path, pos)
% n_chan = number of channels on nirs cap, not number included in analysis
% n_areas = number of Brodmann's areas - eventually we could update this so
% that we can also use other regions of interest, but for now this parameter
% should always be 47
% td_path = file path to find the TDdatabase
% pos = POS file that contains where the probes sat for a single
% participant's session
%% define radius in mm for the sphere around each channel location
rad = 15;
%% define voxel size in mm
vSize = 10;
%% import a template nii and nii resolution
nii = load_nii([nii_path 'avg152T1_LR_nifti.nii']);
vs = nii.hdr.dime.pixdim(2);
%% create the index for voxels
for x = 1:size(nii.img,1)
    for y = 1:size(nii.img,2)
        for z = 1:size(nii.img,3)
            nii.img(x,y,z) = ceil((x-1)*size(nii.img,2)*size(nii.img,3)/vSize) + ...
                ceil((y-1)*size(nii.img,3)/vSize) + ceil(z/vSize);
        end
    end
end
            
%% initialize the transformation matrix
nVoxels = max(nii.img(:));
transformation_mat = zeros(n_chan, nVoxels);

%% then find distance for every channel
channel_locations = pos.R.ch.xyzC';

for chan = 1:n_chan    
    [x,y,z] = mni2xyz(channel_locations(chan,1),channel_locations(chan,2),channel_locations(chan,3));
    chXYZ = [x,y,z];     
    for v = 1:nVoxels
        voxelXYZ = find(nii.img == v);
        transformation_mat(chan, v) = sum(sqrt(sum((voxelXYZ - chXYZ).^2,2)) <= rad/vs);
    end
    transformation_mat(chan, :) = transformation_mat(chan, :)./sum(transformation_mat(chan, :));
end

end



