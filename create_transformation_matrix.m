function mcp_struct = create_transformation_matrix(mcp_struct, convert_to_brodmanns, pathname, probe_loc_file, path_to_tdDatabase, session_idx)
%% create a transformation matrix for each subjects' sessions
% input:
% mcp_struct - mcp data structure for one participant
% convert_to_brodmanns - true or false of whether we want ROI space or
% channel space
% pathname - path to where that participant's session probe location data can be found
%probe_loc_file - name of file where that session's probe locations are
%stored (usually POS.mat)
% path_to_tdDatabase - path to where the TD Database with MNI coordiantes
% for Brodmann's areas can be found

%% create transformation matrices
if convert_to_brodmanns
    % load the POS file
    pos_mat = load([pathname probe_loc_file]);
    n_channels = length(mcp_struct.Experiment.Probe_arrays.Channels);
    transformation_mat = mapChanneltoROI(n_channels, 47, path_to_tdDatabase, pos_mat);
else
    transformation_mat = eye(n_channels);
end

mcp_struct.Experiment.Runs(session_idx).Transformation_Matrix = transformation_mat;


end