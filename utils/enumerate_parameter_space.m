function [parameter_space, user_entered_proc_params, user_entereed_classif_params] = enumerate_parameter_space(results_struct)
%% takes in a struct containing the parameters we want to enumerate in the grid search
% returns a cell array with every combination of those parameters
% first get list of what processing parameters we will search
optional_proc_params = {'scale_data',...
    'scale_function',...
    'incl_sessions',...
    'time_window',...
    'baseline_window',...
    'hemoglobin',...
    'minMax'};

parameter_list = struct;
user_entered_proc_params = {};
user_entereed_classif_params = {};

for p = 1:length(optional_proc_params)
    if isfield(results_struct, optional_proc_params{p}) && iscell(results_struct.(optional_proc_params{p}))
        parameter_list.(optional_proc_params{p}) = results_struct.(optional_proc_params{p});
        user_entered_proc_params{end+1} = optional_proc_params{p};
    end
end

% then get list of what classification parameters we will search
optional_classif_opts = fieldnames(results_struct.opts_struct);
optional_classif_opts_from_user = {};

for p = 1:length(optional_classif_opts)
    if iscell(results_struct.opts_struct.(optional_classif_opts{p}))
        parameter_list.(optional_classif_opts{p}) = results_struct.opts_struct.(optional_classif_opts{p});
        user_entereed_classif_params{end+1} = optional_classif_opts{p};
    end
end


%% now get all combinations of these parameters
fields = fieldnames(parameter_list);
parameter_space = [];

param_fields = fieldnames(parameter_list);
values = struct2cell(parameter_list)';

tmp = cellfun(@(x) 1:length(x), values, 'UniformOutput', false);
combs = combvec(t{:});

params = cell(1,4);

for i = 1:length(values)
    params{i} = values{i}(1,combs(i,:));
end

parameter_space = vertcat(params{:})';

end
















