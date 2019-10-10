function allsubj_results = create_results_struct(summed_mcpa, MCP_dat, parsed_input, mcpa_dat, all_sets, num_subj, num_sets, num_chan, num_cond)

allsubj_results = [];
allsubj_results.patterns = summed_mcpa.patterns;
allsubj_results.MCP_data = MCP_dat;
allsubj_results.created = datestr(now);
allsubj_results.test_handle = parsed_input.Results.test_handle;
allsubj_results.test_type = 'N-fold (Leave one subject out), Classify participant-level averages';
allsubj_results.setsize = parsed_input.Results.setsize;
allsubj_results.func_handle = parsed_input.Results.summary_handle;
allsubj_results.incl_channels = mcpa_dat.incl_channels;
allsubj_results.conditions = parsed_input.Results.conditions;
allsubj_results.subsets = all_sets;

for cond_id = 1:num_cond
    allsubj_results.accuracy(cond_id).condition = allsubj_results.conditions(cond_id);
    allsubj_results.accuracy(cond_id).subjXchan = nan(num_subj,num_chan);
    allsubj_results.accuracy(cond_id).subsetXsubj = nan(num_sets,num_subj);
end

end