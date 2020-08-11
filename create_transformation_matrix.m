function mcp_struct = create_transformation_matrix(mcp_struct, convert_to_ROI, pos_pathname, probe_loc_file, ROIdatabase, areas, session_idx, use_proportional)
%% create a transformation matrix for each subjects' sessions
% input:
% mcp_struct - mcp data structure for one participant
% convert_to_ROI - true or false of whether we want ROI space or
% channel space
% pos_pathname - path to where that participant's session probe location data can be found
% probe_loc_file - name of file where that session's probe locations are
% stored (usually POS.mat)
% ROIdatabase - database with MNI coordiantes where ROIs can be found
% areas - name of ROIs to co-register to as they appear in MNI database
% use_proportional - use proportional assignment instead of assinging one channel to one ROI

%% check that only one subject was entered at a time
if length(mcp_struct) > 1
    error('Only one subject can be entered into create_transformation_matrix at a time. Please specify which participant to use.')
end

%% create transformation matrices
n_channels = length(mcp_struct.Experiment.Probe_arrays.Channels);

if convert_to_ROI
    
    pos_mat = load([pos_pathname probe_loc_file]); % pos file
    transformation_mat = mapChanneltoROI(n_channels, ROIdatabase, pos_mat, areas, use_proportional);

else    
    transformation_mat = eye(n_channels);
end

mcp_struct.Experiment.Runs(session_idx).Transformation_Matrix = transformation_mat;


end
