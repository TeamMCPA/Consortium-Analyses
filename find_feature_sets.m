function sets = find_feature_sets(results)
%% check memory capacity and create all possible subsets
% Arguments:
% results: the results struct from parse_inputs - contains all parameters
    % for classification. Must be passed through as p.Results
    
n_all_sets = nchoosek(length(results.incl_features), results.setsize);
size_of_sets_inmem = n_all_sets*results.setsize*8+100;
%% make sure we don
try
    available_mem = memory;
    available_mem = available_mem.MemAvailableAllArrays;
catch
    try
        [~,w] = unix('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        available_mem = (stats(3)+stats(end))*1000;
    catch
        available_mem = 2E9;
    end
end
%%
if size_of_sets_inmem > 0.50*available_mem
    try
        warning('Too many feature sets will cause memory problems. Randomly generating a subset.');
        sets = nan(results.max_sets,results.setsize);
        for i = 1:results.max_sets
            sets(i,:) = randsample(results.incl_features,results.setsize);
        end
    catch
        error('Too many feature sets will cause memory problems. Reduce setsize or length of incl_features.');
    end
    
else
    if results.verbose
        fprintf('Generating %g possible sets... ',n_all_sets);
    end
    % Start by generating list of all sets
    tic;
    sets = nchoosek(results.incl_features,results.setsize);
    if results.verbose
        fprintf('Done\n');
        toc
    end
    % Case where there are more possible sets than the max_sets limit
    if results.max_sets < n_all_sets
        % Randomly select rows to use
        sets_to_choose = randsample(1:size(sets,1),results.max_sets);
        % Select the randomly sampled rows from sets.
        sets = sets(sets_to_choose,:);
        if results.verbose
            fprintf('Selected %g sets to test.\n',results.max_sets);
        end
    end
end

end






