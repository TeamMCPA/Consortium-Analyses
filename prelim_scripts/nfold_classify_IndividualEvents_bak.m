function allsubj_results = nfold_classify_IndividualEvents(MCP_struct,incl_channels,time_window,cond1,cond2,summary_handle,setsize)

test_handle = @mcpa_trial;

%% Load MCP struct if necessary
if isstring(MCP_struct)
    MCP_struct = load(MCP_struct,'-mat');
    varname = fieldnames(MCP_struct);
    MCP_struct = eval(['MCP_struct.' varname{1}]);
end

%% Build MCPA struct for all subjects in the MCP
% Step 1: Epoching the data by time window and averaging the epochs
% together at the subject level
mcpa_struct = MCP_to_MCPA(MCP_struct,[1:length(MCP_struct)],incl_channels,time_window);
% Step 2: Apply the desired function (e.g., @nanmean) for summarizing time
% window data. You can write custom functions to deal with time- and
% channel-domain data however you want. Default behavior is to apply the
% function along the first dimension of the MCPA pattern, but this can also
% be changed.
mcpa_summ = summarize_MCPA_Struct(summary_handle,mcpa_struct);

%% Identify the conditions that will be compared
% Condition flags can be string, cell, integer, or logical. If string or
% cell, a search is run over the summarized MCPA data to determine which
% condition # matches the condition name. Otherwise, integer or logical
% values are applied directly to indexing.
if ischar(cond1) || isstring(cond1) || iscellstr(cond1)
    cond1_flag = strcmp(cond1,mcpa_summ.event_types);
else
    cond1_flag = cond1;
end
if ischar(cond2) || isstring(cond2) || iscellstr(cond2)
    cond2_flag = strcmp(cond2,mcpa_summ.event_types);
else
    cond2_flag = cond2;
end

%% Setting up the combinations of channel subsets
% Create all possible subsets. If setsize is equal to the total number of
% channels, there will only be one 'subset' which is the full channel
% array. If setsize is less than the total number of channels, there will
% be n-choose-k subsets to analyze.
sets = nchoosek(incl_channels,setsize);

%% Set up the results structure which includes a copy of MCPA_pattern
allsubj_results = [];
allsubj_results.MCPA_patterns = mcpa_struct.patterns;
allsubj_results.MCP_data = MCP_struct;
allsubj_results.created = datestr(now);
allsubj_results.test_handle = test_handle;
allsubj_results.test_type = 'N-fold (Leave one subject out), Classify individual events';
allsubj_results.setsize = setsize;
allsubj_results.func_handle = summary_handle;
allsubj_results.incl_channels = mcpa_struct.incl_channels;
allsubj_results.cond1 = cond1;
allsubj_results.cond2 = cond2;
allsubj_results.subsets = sets;

n_subj = length(MCP_struct);
n_sets = size(sets,1);
n_chan = length(incl_channels);
n_events = max(arrayfun(@(x) max(sum(x.fNIRS_Data.Onsets_Matrix)),MCP_struct));

allsubj_results.accuracy.cond1.subjXchan = nan(n_subj,n_chan);
allsubj_results.accuracy.cond2.subjXchan = nan(n_subj,n_chan);
allsubj_results.accuracy.cond1.subsetXsubj = nan(n_sets,n_subj);
allsubj_results.accuracy.cond2.subsetXsubj = nan(n_sets,n_subj);

for s_idx = 1:length(MCP_struct)
    
    fprintf('Running %g feature subsets for Subject %g / %g',n_sets,s_idx,n_subj);
    
    %% Extract training data
    group_subvec = 1:length(MCP_struct);
    group_subvec(s_idx) = [];
    
    group_cond1 = mean(mcpa_summ.patterns(cond1_flag,incl_channels,group_subvec),3)';
    group_cond2 = mean(mcpa_summ.patterns(cond2_flag,incl_channels,group_subvec),3)';
    
    
    %% Extract test data
    subj_events = MCP_get_subject_events(MCP_struct(s_idx),incl_channels,time_window,{MCP_struct(s_idx).Experiment.Conditions.Name});
    subj_patterns1 = squeeze(summary_handle(subj_events(:,:,:,cond1_flag)));
    subj_patterns2 = squeeze(summary_handle(subj_events(:,:,:,cond2_flag)));
    

    %% Run over channel subsets
    tic;
    for set_idx = 1:n_sets
        
        % Report at every 5% progress
        status_jump = floor(n_sets/20);
        if ~mod(set_idx,status_jump)
            fprintf(' .')
        end
        
        temp_acc1 = nan(n_events,1);
        temp_acc2 = nan(n_events,1);
        
        set_chans = sets(set_idx,:);
        for e_idx = 1:size(subj_patterns1,2)
            temp_acc1(e_idx) = corr(subj_patterns1(set_chans,e_idx),group_cond1(set_chans)) > corr(subj_patterns1(set_chans,e_idx),group_cond2(set_chans));
            %temp_acc1(e_idx) = ( test_handle(subj_patterns1(set_chans,e_idx), group_cond1(set_chans), group_cond2(set_chans))==1 );
        end
        for e_idx = 1:size(subj_patterns2,2)
            temp_acc2(e_idx) = corr(subj_patterns2(set_chans,e_idx),group_cond2(set_chans)) > corr(subj_patterns2(set_chans,e_idx),group_cond1(set_chans));
            %temp_acc2(e_idx) = ( test_handle(subj_patterns1(set_chans,e_idx), group_cond2(set_chans), group_cond1(set_chans))==2 );
        end
        
        allsubj_results.accuracy.cond1.subjXchan(s_idx,set_chans) = nanmean(temp_acc1);
        allsubj_results.accuracy.cond2.subjXchan(s_idx,set_chans) = nanmean(temp_acc2);
        allsubj_results.accuracy.cond1.subsetXsubj(set_idx,s_idx) = nanmean(temp_acc1);
        allsubj_results.accuracy.cond2.subsetXsubj(set_idx,s_idx) = nanmean(temp_acc2);
    end
    fprintf('%0.1f mins\n',toc/60);
    
end

function classification = mcpa_trial(obsv_pattern, model_pattern1, model_pattern2)
if corr(obsv_pattern,model_pattern1) > corr(obsv_pattern,model_pattern2)
    classification = 1;
else
    classification = 2;
end
end

%% Visualization
if n_sets > 1

    figure
    errorbar(1:size(allsubj_results.accuracy.cond1.subjXchan,2),mean(allsubj_results.accuracy.cond1.subjXchan),std(allsubj_results.accuracy.cond1.subjXchan)/sqrt(size(allsubj_results.accuracy.cond1.subjXchan,1)),'r')
    hold;
    errorbar(1:size(allsubj_results.accuracy.cond2.subjXchan,2),mean(allsubj_results.accuracy.cond2.subjXchan),std(allsubj_results.accuracy.cond2.subjXchan)/sqrt(size(allsubj_results.accuracy.cond2.subjXchan,1)),'k')
    title('Decoding Accuracy across all channels: Red = Cond1, Black = Cond2')
    set(gca,'XTick',[1:length(incl_channels)])
    set(gca,'XTickLabel',incl_channels)
    hold off;
    
    figure
    errorbar(1:size(allsubj_results.accuracy.cond1.subjXchan,1),mean(allsubj_results.accuracy.cond1.subjXchan'),repmat(std(mean(allsubj_results.accuracy.cond1.subjXchan'))/sqrt(size(allsubj_results.accuracy.cond1.subjXchan,2)),1,size(allsubj_results.accuracy.cond1.subjXchan,1)),'r')
    hold;
    errorbar(1:size(allsubj_results.accuracy.cond2.subjXchan,1),mean(allsubj_results.accuracy.cond2.subjXchan'),repmat(std(mean(allsubj_results.accuracy.cond2.subjXchan'))/sqrt(size(allsubj_results.accuracy.cond2.subjXchan,2)),1,size(allsubj_results.accuracy.cond2.subjXchan,1)),'k')
    title('Decoding Accuracy across all subjects: Red = Cond1, Black = Cond2')
    set(gca,'XTick',[1:length(MCP_struct)])
    hold off;
    
end

end