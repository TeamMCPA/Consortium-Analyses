function [scaled_full_set, scaled_testing_set, scaled_training_set] = standardize_data(full_dataset, testing_data, training_data, input_struct)
%% standardize data so that we are working with z-scores
% Arguments:
% group: group model data (training data for classification)
% subject: left out subject data (testing data for classification)

scaled_full_set = [];
scaled_testing_set = [];
scaled_training_set = [];


if input_struct.norm_features % feature-wise scaling
    if isempty(full_dataset) % scale split data
        average = mean(training_data,1);
        sd = std(training_data,0,1);
        scaled_training_set = (training_data - average) ./ sd;
        scaled_testing_set = (testing_data - average)./ sd;
    else % scale individual data
        average = mean(full_dataset,1);
        sd = std(full_dataset,0,1);
        scaled_full_set = (full_dataset - average) ./ sd;
    end  
else % row-wise scaling
    if isempty(full_dataset) % scale split data
        average = mean(training_data,2);
        sd = std(training_data,0,2);
        scaled_training_set = (training_data - average) ./ sd;
        scaled_testing_set = (testing_data - average)./ sd;
    else % scale individual data
        average = mean(full_dataset,2);
        sd = std(full_dataset,0,2);
        scaled_full_set = (full_dataset - average) ./ sd;
    end  
end

end