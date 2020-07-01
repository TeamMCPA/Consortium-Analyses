function mcp_struct = create_transformation_matrix(mcp_struct, convert_to_brodmanns, split_hemispheres, pos_pathname, probe_loc_file, tdDatabase, session_idx)
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

n_channels = length(mcp_struct.Experiment.Probe_arrays.Channels);

if convert_to_brodmanns
    
    pos_mat = load([pos_pathname probe_loc_file]); % pos file
    areas = {};
    for area = 1:47
        area_name = ['brodmann_area_', int2str(area)];
        areas{end+1} = area_name;
    end
    
    if split_hemispheres
        brodmanns_areas_hemispheres = struct;
        
        for a = 1:length(areas)
            temp_area = tdDatabase.(areas{a});
            left_area = find(temp_area(:,1) < 0);
            right_area = find(temp_area(:,1) > 0);
            new_name_left = [areas{a} '_left'];
            new_name_right = [areas{a} '_right'];
            brodmanns_areas_hemispheres.(new_name_left) = temp_area(left_area,:);
            brodmanns_areas_hemispheres.(new_name_right) = temp_area(right_area,:);
        end
        
        areas = fieldnames(brodmanns_areas_hemispheres);
        
        transformation_mat = mapChanneltoROI(n_channels, brodmanns_areas_hemispheres, pos_mat, areas);
        
    else
        transformation_mat = mapChanneltoROI(n_channels, tdDatabase, pos_mat,areas);
   
    end
else
    
    transformation_mat = eye(n_channels);
end

mcp_struct.Experiment.Runs(session_idx).Transformation_Matrix = transformation_mat;


end
