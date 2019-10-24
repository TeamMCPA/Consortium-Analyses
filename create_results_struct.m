function allsubj_results = create_results_struct(summed_mcpa, MCP_dat, parsed_input, mcpa_dat, all_sets, num_subj, num_sets, num_chan, num_cond)
%% creates a results struct to output classification accuracies to 
% Arguments:
% summed_mcpa: summarized MCPA struct from Summarize_MCPA_Struct
% MCP_dat: MCP struct
% parsed_input: struct containing either default classification parameters 
% or those passed through when calling nfold_classify
% all_sets: all channel subsets
% num_subj: number of subjects
% num_sets: number of sets we use for classification 
% num_chan: number of channels
% num_cond: number of conditions


allsubj_results = []; % create empty structs
allsubj_results.patterns = summed_mcpa.patterns; % store the multichannel patterns from the summarized MCPA
allsubj_results.MCP_data = MCP_dat; % store MCP struct
allsubj_results.created = datestr(now); % date the struct is created
allsubj_results.test_handle = parsed_input.Results.test_handle; % what method of classification did we use Default: mcpa_classify
allsubj_results.test_type = 'N-fold (Leave one subject out), Classify participant-level averages'; % what type of test we used
allsubj_results.setsize = parsed_input.Results.setsize; % size of channel sets
allsubj_results.func_handle = parsed_input.Results.summary_handle; % what method did we use to summarize our data Default: @nanmean
allsubj_results.incl_channels = mcpa_dat.incl_channels; % what channels did we include
allsubj_results.conditions = parsed_input.Results.conditions; % list of conditions
allsubj_results.subsets = all_sets; % all channel subsets

for cond_id = 1:num_cond % now create place holders for decoding accuracies 
    allsubj_results.accuracy(cond_id).condition = allsubj_results.conditions(cond_id);
    allsubj_results.accuracy(cond_id).subjXchan = nan(num_subj,num_chan);
    allsubj_results.accuracy(cond_id).subsetXsubj = nan(num_sets,num_subj);
end

end