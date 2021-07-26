function semantic_model = validate_model_matrix(semantic_model, input_struct)
%% check that the model matrix is in fact already in representational (dis)similarity space
% if not, abstract it into representational (dis)similarity space using the
% same metrics as for the held out subject's data

if size(semantic_model,1) ~= size(semantic_model,2)    
    semantic_model = input_struct.opts.similarity_function(model_dat, input_struct.opts);
end



end