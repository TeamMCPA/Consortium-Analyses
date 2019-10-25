function [group_data, group_labels, subj_data, subj_labels] = split_data(s_idx, cond_flags, p, n_cond, mcpa_summ,group_subvec, group_data, group_labels, subj_data, subj_labels)	
%% Split the training and testing data and create vectors for their labels
% Arguments:
% s_idx: which subject we left out
% cond_flags: flags for what condition we are using
% p: input struct from parse_inputs
% n_cond: number of conditions we are using
% mcpa_summ: summarized MCPA struct (multichannel patterns)
% group_subvec: indices for which subjects are included in the group data
% group_data: training data
% group_labels: training labels
% subj_data: testing data
% subj_labels: testing labels

for cond_idx = 1:n_cond
        if ischar(p.Results.conditions{cond_idx}) || isstring(p.Results.conditions{cond_idx}) || iscellstr(p.Results.conditions{cond_idx})
            [~, ~, cond_flags{cond_idx}] = intersect(p.Results.conditions{cond_idx},mcpa_summ.event_types);
        else
            cond_flags{cond_idx} = p.Results.conditions{cond_idx};
        end
            
	% Extract training data
	% group_data_tmp averages across all matching triggers for a
	% condition and outputs a subj-x-chan matrix
	group_data_tmp = squeeze(mean(mcpa_summ.patterns(cond_flags{cond_idx},p.Results.incl_channels,group_subvec),1))';
	group_labels_tmp = repmat(cellstr(strjoin(string(p.Results.conditions{cond_idx}),'+')),size(group_data_tmp,1),1);
	group_data = [ group_data; group_data_tmp ];
	group_labels = [ group_labels; group_labels_tmp ];
            
	% Extract test data
	subj_data_tmp = mcpa_summ.patterns(cond_flags{cond_idx},p.Results.incl_channels,s_idx);
	subj_labels_tmp = repmat(cellstr(strjoin(string(p.Results.conditions{cond_idx}),'+')),size(subj_data_tmp,1),1);
	subj_data = [ subj_data; subj_data_tmp ];
	subj_labels = [ subj_labels; subj_labels_tmp ];
    end
end
