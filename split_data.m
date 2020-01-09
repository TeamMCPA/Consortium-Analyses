function [group_data, group_labels, subj_data, subj_labels] = split_data(fold_idx, cond_flags, Results, num_cond, patterns,group_subvec, group_data, group_labels, subj_data, subj_labels, event_types)	

if length(size(patterns)) > 3
    length(size(patterns))
    pattern_dat = create_instanceXfold_dimension(patterns);
    repetitions = size(patterns,3);
    start_remove_idx = (repetitions * fold_idx) - repetitions + 1; 
    stop_remove_idx = repetitions * fold_idx;
    
    group_subvec = [1:size(pattern_dat,3)];
    group_subvec(start_remove_idx:stop_remove_idx) = [];
    
    fold_idx = start_remove_idx:stop_remove_idx;    
else
    pattern_dat = patterns;
end


% first split up what data needs to go into testing and training
for cond_idx = 1:num_cond
    if ischar(Results.conditions{cond_idx}) || isstring(Results.conditions{cond_idx}) || iscellstr(Results.conditions{cond_idx})
        [~, ~, cond_flags{cond_idx}] = intersect(Results.conditions{cond_idx},event_types);
    else
        cond_flags{cond_idx} = Results.conditions{cond_idx};
    end
    
    
    % Extract training data
	% group_data_tmp averages across all matching triggers for a
	% condition and outputs a subj-x-chan matrix
    group_data_tmp = squeeze(mean(pattern_dat(cond_flags{cond_idx},Results.incl_channels,group_subvec),1))';
	group_labels_tmp = repmat(cellstr(strjoin(string(Results.conditions{cond_idx}),'+')),size(group_data_tmp,1),1);
	group_data = [ group_data; group_data_tmp ];
	group_labels = [ group_labels; group_labels_tmp ];

    % Extract test data
	subj_data_tmp = squeeze(mean(pattern_dat(cond_flags{cond_idx},Results.incl_channels,fold_idx),1))';
	subj_labels_tmp = repmat(cellstr(strjoin(string(Results.conditions{cond_idx}),'+')),size(subj_data_tmp,1),1);
	subj_data = [ subj_data; subj_data_tmp ];
	subj_labels = [ subj_labels; subj_labels_tmp ];

end
    
    
end
