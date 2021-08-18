function [semantic_model,semantic_model_labels] = validate_model_matrix(semantic_model, semantic_model_labels, conditions, input_struct)
%% check that the model matrix is in fact already in representational (dis)similarity space
% if not, abstract it into representational (dis)similarity space using the
% same metrics as for the held out subject's data

if size(semantic_model,1) ~= size(semantic_model,2)    
    semantic_model = input_struct.opts.similarity_function(model_dat, input_struct.opts);
end

%% then make sure only the needed conditions are included

if size(semantic_model,1) ~= length(conditions)
    % Set logical flags for indexing the conditions that will be compared.
    % Loop through the whole list of conditions and create flags for each.
    cond_flags = cell(length(conditions),1);
    
    for cond_idx = 1:length(conditions)
        if ischar(conditions{cond_idx}) || isstring(conditions{cond_idx}) || iscellstr(conditions{cond_idx})
            [~, ~, cond_flags{cond_idx}] = intersect(conditions{cond_idx},semantic_model_labels);
        else
            cond_flags{cond_idx} = conditions{cond_idx};
        end
    end
    cond_flags = [cond_flags{:}]';
    
    semantic_model = semantic_model(cond_flags, cond_flags);
    semantic_model_labels = semantic_model_labels(cond_flags);
    
end


end
