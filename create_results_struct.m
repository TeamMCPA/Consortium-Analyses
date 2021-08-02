function allsubj_results = create_results_struct(parsed_input, cv_function, all_sets, num_subj, num_sets, num_feature, num_cond, current_dimensions, final_dimensions, pattern_data, event_types)
%% create a struct to store classification results. 
% This will be later used in permutation testing 

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
allsubj_results.event_types = event_types;
allsubj_results.test_type = cv_function;


%% create accuracy struct
for cond_id = 1:num_cond % now create place holders for decoding accuracies 
    allsubj_results.accuracy(cond_id).condition = allsubj_results.conditions(cond_id);
    allsubj_results.accuracy(cond_id).subjXfeature = nan(num_subj,num_feature);
    allsubj_results.accuracy(cond_id).subsetXsubj = nan(num_sets,num_subj);
    if strcmp(cv_function, 'nfold_classify_WithinSubjects')
         allsubj_results.accuracy(cond_id).subjXsession = nan(num_subj,4);
    end
end

if strcmp(cv_function, 'nfold_generalize_ParticipantLevel')
    % This renames the conditions for the accuracy field - currently create_results_struct operates as
    % though we're using all the conditions so it labels the accuracy fields as
    % baby 1, baby 2, etc. when we really want baby, bottle, etc.
    groups = unique(parsed_input.cond_key(:,2));
    for group_id = 1:length(groups)
        allsubj_results.accuracy(group_id).condition = groups(group_id);
    end
end

%% create an accuracy matrix if doing pairwise
if strcmp(cv_function, 'nfold_classify_WithinSubjects')
    if isfield(parsed_input.opts_struct, 'pairwise') && parsed_input.opts_struct.pairwise == 1
        allsubj_results.accuracy_matrix = nan(length(allsubj_results.conditions),...
            length(allsubj_results.conditions),...
            min(size(allsubj_results.subsets,1),...
            allsubj_results.max_sets),...
            4,...
            num_subj);
    end
else
    if isfield(parsed_input.opts_struct, 'pairwise') && parsed_input.opts_struct.pairwise == 1
        allsubj_results.accuracy_matrix = nan(length(allsubj_results.conditions),...
            length(allsubj_results.conditions),...
            min(size(allsubj_results.subsets,1),...
            allsubj_results.max_sets),...
            num_subj);
    end

end
