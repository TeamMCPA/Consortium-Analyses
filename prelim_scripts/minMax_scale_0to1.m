function [group_dat, test_dat] = normalize_data(group, test)

mins = min(group);
maxes = max(group);
group_dat = (group - mins) ./ (maxes - mins);
test_dat = (test-mins) ./ (maxes-mins);


end