function temp_results_struct = populate_classification_options_struct(results_struct, parameter_space, proc_params, classif_params, p_idx, is_inner_loop)

temp_results_struct = results_struct;
for p = 1:length(proc_params)
    temp_results_struct.(proc_params{p}) = parameter_space{p_idx, p};
end

opts_struct = results_struct.opts_struct;
for p = 1:length(classif_params)
    opts_struct.(classif_params{p}) = parameter_space{p_idx, p+length(proc_params)};
end


temp_results_struct.opts_struct = opts_struct;

if is_inner_loop
    temp_results_struct.incl_subjects = 1:(length(temp_results_struct.incl_subjects)-1);
end

end


