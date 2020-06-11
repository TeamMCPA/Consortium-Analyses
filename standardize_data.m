function [group_dat, subject_dat] = standardize_data(group, subject)
%% standardize data so that we are working with z-scores
% Arguments:
% group: group model data (training data for classification)
% subject: left out subject data (testing data for classification)

average = mean(group);
sd = std(group);
group_dat = (group-average) ./ sd;

% only want to standardize the subject data if we're trying to do this within cross validation. 
% If we're doing this in CV,we want to apply the group's sd and mean to the subject's data
if ~isstruct(subject)
    subject_dat = (subject - average) ./ sd;
else
    subject_dat = nan;
end

end
