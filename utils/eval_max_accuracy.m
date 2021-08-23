function [chosen_parameters, rankings] = eval_max_accuracy(grid_search_accuracies, parameter_space)
%% evaluate grid search results by finding the maximum accuracy
% input: 
% grid_search_accuracies - vector of the accuracies in each run of the grid
% search
% parameter_space - every combination of parameters
% returns the parameters that lead to the highest accuracy and a ranking of
% how every other parameter performed

parameter_combinations = fieldnames(grid_search_accuracies);

acc_from_parameterset = [];
for param = 1:size(parameter_combinations,1)
    if isfield(grid_search_accuracies.(parameter_combinations{param}),'accuracy_matrix')
        acc_matrix = grid_search_accuracies.(parameter_combinations{param}).accuracy_matrix;
        participant_accs = [];
        for i = 1: (size(acc_matrix,3)*size(acc_matrix,4))            
            tmp = acc_matrix(:,:,i);
            participant_accs = [participant_accs; nanmean(tmp(:))];
        end
        acc_from_parameterset = [acc_from_parameterset; nanmean(participant_accs)];
            
    else    
        acc_from_conditions = [];
        for cond = 1:length(grid_search_accuracies.(['Parameter_set_' int2str(param)]).accuracy)
            acc_from_conditions = [acc_from_conditions; nanmean(nanmean(grid_search_accuracies.(['Parameter_set_' int2str(param)]).accuracy(cond).subsetXsubj))];
        end

        acc_from_parameterset = [acc_from_parameterset; nanmean(acc_from_conditions)];
    end
    
end


[~, rankings] = sort(acc_from_parameterset,'descend');
chosen_parameters = parameter_space(rankings(1),:);

end