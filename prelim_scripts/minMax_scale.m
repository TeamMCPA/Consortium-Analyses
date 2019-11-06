function [scaled_full_set, scaled_testing_set, scaled_training_set] = minMax_scale(full_dataset, testing_data, training_data, input_struct)
%% This function scales data to a set minimum and maximum
% Arguments: 
% full_dataset: Either input a full Oxy dataset from one participant if you
% want to scale individual subject data or enter this as [] if you want to
% scale group level data
% testing_data: Either input the testing data (after we split test and
% train data) or input [] if we are scaling individual level data
% training_data: Input training data (after we split test and train data)
% or input [] if scaling individual level data
% input_struct: stucture we entered into folding function with our
% parameters


minimum = input_struct.minMax(1); % what minimum are we scaling to
maximum = input_struct.minMax(2); % what maximum are we scaling to

scaled_full_set = [];
scaled_testing_set = [];
scaled_training_set = [];

if input_struct.norm_features % feature-wise scaling
    if isempty(full_dataset) % scale split data
        mins = min(training_data);
        maxes = max(training_data);
        scaled_training_set = (maximum - minimum)*((training_data - mins) ./ (maxes - mins)) + minimum;
        scaled_testing_set = (testing_data-mins) ./ (maxes-mins);
    else % scale individual data
        mins = min(full_dataset);
        maxes = max(full_dataset);
        scaled_full_set = (maximum - minimum)*((full_dataset - mins) ./ (maxes - mins)) + minimum;
    end
else % row-wise scaling
    if isempty(full_dataset)
        mins = min(training_data, [], 2);
        maxes = max(training_data, [], 2);
        scaled_training_set = (maximum - minimum)*((training_data - mins) ./ (maxes - mins)) + minimum;
        scaled_testing_set = (testing_data-mins) ./ (maxes-mins);
    else
        mins = min(full_dataset, [], 2);
        maxes = max(full_dataset, [], 2);
        scaled_full_set = (maximum - minimum)*((full_dataset - mins) ./ (maxes - mins)) + minimum;
    end

end
end