function [group_dat, subject_dat] = minMax_scale_0to1(group, subject)
%% resacle the data so that the minimum is 0 and the maximum is 1
% Arguments:
% group: group model data (training data for classification)
% subject: left out subject data (testing data for classification)

mins = min(group);
maxes = max(group);
group_dat = (group - mins) ./ (maxes - mins);
subject_dat = (subject-mins) ./ (maxes-mins);


end