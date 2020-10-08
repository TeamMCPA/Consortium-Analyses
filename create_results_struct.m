function allsubj_results = create_results_struct(parsed_input, cv_function, all_sets, num_subj, num_sets, num_feature, num_cond, current_dimensions, final_dimensions, pattern_data)
%% create a struct to store classification results. 
% This will be later used in permutation testing, cross validation, and nested cross validation

% created by Anna Herbolzheimer and Ben Zinszer 2019
% revised by Anna Herbolzheimer fall 2020

allsubj_results = struct; % create empty structs
allsubj_results.created = datestr(now);

input_fields = fields(parsed_input);
input_values = struct2cell(parsed_input);

% fill the results struct with our parameters
for f = 1:length(input_fields)
    allsubj_results.(input_fields{f}) = input_values{f};
end

allsubj_results.final_dimensions = final_dimensions;
allsubj_results.subsets = all_sets;
allsubj_results.patterns = pattern_data;
allsubj_results.dimensions = current_dimensions;



%% create accuracy struct
for cond_id = 1:num_cond % now create place holders for decoding accuracies 
    allsubj_results.accuracy(cond_id).condition = allsubj_results.conditions(cond_id);
    allsubj_results.accuracy(cond_id).subjXfeature = nan(num_subj,num_feature);
    allsubj_results.accuracy(cond_id).subsetXsubj = nan(num_sets,num_subj);
    if strcmp(func2str(cv_function), 'cross_validate_WithinSubjects')
         allsubj_results.accuracy(cond_id).subjXsession = nan(num_subj,max_sessions);
    end
end

%% create an accuracy matrix if doing pairwise classification
if isfield(parsed_input.opts_struct.pairwise) && parsed_input.opts_struct.pairwise == 1
    allsubj_results.accuracy_matrix = nan(length(allsubj_results.conditions),...
        length(allsubj_results.conditions),...
        min(size(allsubj_results.subsets,1),...
        allsubj_results.max_sets),...
        size(pattern_data, ndims(pattern_data)-1),...
        size(pattern_data, ndims(pattern_data)));
end
