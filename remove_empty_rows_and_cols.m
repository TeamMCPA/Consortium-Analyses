function [model_dat, test_dat] = remove_empty_rows_and_cols(model_data, test_data, model_dat, test_dat)
%% first for the model data
empty_x_vals_model = [];
empty_y_vals_model = [];
for i = 1:size(model_data,4)
    for j = 1:size(model_data,3)
        [x,y] = find(isnan(model_dat(:,:,j,i)));
        if length(unique(x)) == size(model_dat,1) && length(unique(y)) == size(model_dat,2)
            continue;
        else
            empty_x_vals_model = [empty_x_vals_model; x];
            empty_y_vals_model = [empty_y_vals_model; y];
        end
    end      
end
empty_x_vals_model = unique(empty_x_vals_model);
empty_y_vals_model = unique(empty_y_vals_model);

% then test data
empty_x_vals_test = [];
empty_y_vals_test = [];
for i = 1:size(test_data,4)
    for j = 1:size(test_data,3)
        [x,y] = find(isnan(test_data(:,:,j,i)));
        
        
        if length(unique(x)) == size(test_data,1) && length(unique(y)) == size(test_data,2)
            continue;
        else
            empty_x_vals_test = [empty_x_vals_test; x];
            empty_y_vals_test = [empty_y_vals_test; y];
        end
    end      
end


empty_x_vals_test = unique(empty_x_vals_test);
empty_y_vals_test = unique(empty_y_vals_test);

cols = 1:size(model_data,2);
rows = 1:size(model_data,1);

if size(model_data,2) > 8
    if length(empty_x_vals_model) == size(model_data,1) && length(empty_y_vals_model) ~= size(model_data,2)
        %remove columns
        remove_cols = union(empty_y_vals_model, empty_y_vals_test);
        keep = ~ismember(cols, remove_cols);
        model_dat = model_dat(:,keep', :,:);
        test_dat = test_dat(:, keep', :,:);   
    elseif length(empty_x_vals_model) == size(model_data,1) && length(empty_y_vals_model) == size(model_data,2)
        %this is an empty session
        warning('There is no data in this session')
    end
else
    if length(empty_x_vals_test) == size(model_data,1) && length(empty_y_vals_test) ~= size(test_data,2)
        test_cols = 1:size(test_dat,2);
        remove_cols = union(empty_y_vals_model, empty_y_vals_test);
        keep = ~ismember(test_cols, remove_cols);
        test_dat = test_dat(:, keep', :,:);
    end 
end

end