function [simulated_data,simulated_labels] = simulate(channel_category_means, channel_category_sds, time_sd, num_chans, num_trials, num_categories)
%% function to simulate NIRS data
% Input:
% channel_category_means - a vector containing the mean value for each
% channel, e.g. [10, 30] or [10, 10]
% channel_category_sds - a vector containing the standard deviation for
% each channel 
% num_chans - number of channels
% num_trials - number of trials for each category - if classes are
% balanced, input 1 value, otherwise enter a vector of how many trials
% exist for each class
% num_categories - number of categories

temp_simulated_data = nan(max(num_trials), num_chans,num_categories);
simulated_labels = [];

for cat = 1:num_categories
    
    if length(num_trials) ~= 1
        n_trials_cat = num_trials(cat);
    else
        n_trials_cat = num_trials;
    end
    
    ch_mean = normrnd(channel_category_means(cat), channel_category_sds(cat), [1,num_chans]);
    
    for chan = 1:num_chans
        ch_across_time = normrnd(ch_mean(chan),time_sd, [1,n_trials_cat]);
        temp_simulated_data(1:n_trials_cat, chan, cat) = ch_across_time';
    end
    
    simulated_labels = [simulated_labels; repmat(cat, n_trials_cat, 1)];
    
end


simulated_data=concatenate_dimensions(temp_simulated_data, [1,3]);





end
