function [group_dat, test_dat] = standardize_data(group, test)

average = mean(group);
sd = std(group);
group_dat = (group-average) ./ sd;
test_dat = (test - average) ./ sd;

end