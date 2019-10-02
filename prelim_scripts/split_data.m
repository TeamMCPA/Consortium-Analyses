function [group_dat, group_labs, test_dat, test_labs] = split_data(group_dat, group_labs, test_dat, test_labs, num_cond, condition_flags, parsed_input, summed_mcpa, g_subvec, sub_idx)

for cond_idx = 1:num_cond
    if ischar(parsed_input.Results.conditions{cond_idx}) || isstring(parsed_input.Results.conditions{cond_idx}) || iscellstr(parsed_input.Results.conditions{cond_idx})
        [~, ~, condition_flags{cond_idx}] = intersect(parsed_input.Results.conditions{cond_idx},summed_mcpa.event_types);
    else
        condition_flags{cond_idx} = parsed_input.Results.conditions{cond_idx};
    end
            
    % Extract training data
    % group_data_tmp averages across all matching triggers for a
    % condition and outputs a subj-x-chan matrix
	group_data_tmp = squeeze(mean(summed_mcpa.patterns(condition_flags{cond_idx},parsed_input.Results.incl_channels,g_subvec),1))';
	group_labels_tmp = repmat(cellstr(strjoin(string(parsed_input.Results.conditions{cond_idx}),'+')),size(group_data_tmp,1),1);
	group_dat = [ group_dat; group_data_tmp ];
	group_labs = [ group_labs; group_labels_tmp ];
            
	% Extract test data
	subj_data_tmp = summed_mcpa.patterns(condition_flags{cond_idx},parsed_input.Results.incl_channels,sub_idx);
	subj_labels_tmp = repmat(cellstr(strjoin(string(parsed_input.Results.conditions{cond_idx}),'+')),size(subj_data_tmp,1),1);
	test_dat = [ test_dat; subj_data_tmp ];
	test_labs = [ test_labs; subj_labels_tmp ];
end
end



