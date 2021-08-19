function [simulated_data,simulated_labels] = simulate_data(channel_category_means, channel_category_sds, time_sds, num_chans, num_trials, num_categories)
%% function to simulate NIRS data
% Input:
% channel_category_means - a vector containing the mean value for each
% channel, e.g. [10, 30] or [10, 10]
% channel_category_sds - a vector containing the standard deviation for
% each channel 
% time_sds - a vector containing the standard deviation over time for each category 
% num_chans - number of channels
% num_trials - number of trials for each category - if classes are
% balanced, input 1 value. Otherwise enter a vector of how many trials
% exist for each category
% num_categories - number of categories

simulated_data = nan(max(num_trials), num_chans,num_categories);
simulated_labels = [];

for cat = 1:num_categories
    
    if length(num_trials) ~= 1 % check to see if user input a different amount of trials for each category
        n_trials_cat = num_trials(cat); % if so, increment through the number of trials
    else
        n_trials_cat = num_trials; % otherwise don't increment through the number of trials
    end
    
    % simulate the mean of each channel for that category
    chan_means = normrnd(channel_category_means(cat), channel_category_sds(cat), [1,num_chans]);
    
    for chan = 1:num_chans
        % simulate full time series for each channel
        chan_across_time = normrnd(chan_means(chan),time_sds(cat), [1,n_trials_cat]);
        simulated_data(1:n_trials_cat, chan, cat) = chan_across_time';
    end
    
    simulated_labels = [simulated_labels; repmat(cat, n_trials_cat, 1)];
    
end


simulated_data = concatenate_dimensions(simulated_data, [1,3]); % concatenate_dimensions is a function inside of our decoding toolbox

% after reshape, alternate trials for categories 1, 2, etc. 
key = repmat([1:num_trials]', num_categories,1);
[~, order] = sort(key);
simulated_data=simulated_data(order,:);
simulated_labels=int2str(simulated_labels(order,:));


end
