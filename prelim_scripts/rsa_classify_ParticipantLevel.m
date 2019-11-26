function allsubj_results = rsa_classify_ParticipantLevel(MCP_struct, semantic_model, varargin)
%% rsa_classify_IndividualSubjects takes an MCP struct and performs
% RSA classification for n subjects to classify individual
% participants' average response patterns using a semantic model. This wrapper assumes that
% features will be averaged within-participants to produce a single
% participant-level observation. Thus the training set is constrained to
% the number of participants minus 1. Several parameters can be changed,
% including which functions are used to generate features and what
% classifier is trained. See Arguments below:
%
% Arguments:
% MCP_struct: either an MCP-formatted struct or the path to a Matlab file
% (.mat or .mcp) containing the MCP_struct.
% semantic_model: The model we want to use for classification
% incl_channels: channels to include in the analysis. Default: all channels
% incl_subjects: index of participants to include. Default: all participants
% baseline_window: [onset, offset] in seconds. Default [-5,0]
% time_window: [onset, offset] in seconds. Default [2,6]
% conditions: cell array of condition names / trigger #s. Default: {1,2}
% summary_handle: function handle (or char of function name) to specify how
% time-x-channel data should be summarized into features. Default: nanmean
% setsize: number of channels to analyze (for subset analyses) Default: all
% test_handle: function handle for classifier. Default: mcpa_classify
% opts_struct: contains additional classifier options. Default: empty struct
% verbose: logical flag to report status updates and results. Default: true

%% Load MCP struct if necessary
if isstring(MCP_struct) || ischar(MCP_struct)
    MCP_struct = load(MCP_struct,'-mat');
    varname = fieldnames(MCP_struct);
    MCP_struct = eval(['MCP_struct.' varname{1}]);
end

%% Parse out the input data
p = parse_inputs(MCP_struct, varargin{:});

%% Setting up the combinations of channel subsets
% Create all possible subsets. If setsize is equal to the total number of
% channels, there will only be one 'subset' which is the full channel
% array. If setsize is less than the total number of channels, there will
% be n-choose-k subsets to analyze.
%
% The size of the subsets can grow extremely quickly with the size of
% incl_channels. Consequently, there is a default max of 1000000 sets,
% which can be customized. If the total number of sets is larger than the
% max number of allowed sets, the list of sets will be subsampled.

% Determine how many sets will be generated. Can use this later for warning
% messages or other branching. Sets variable turns into a huge memory hog.
unmapped_sets = find_sets(p.Results);
sets = map_values(p, unmapped_sets);

%% norm check - do we want to scale individual participant data?
if p.Results.norm_data_participantLevel
    MCP_struct = scale_individuals(MCP_struct, p.Results, sets);
end
%% Build MCPA struct for all subjects in the MCP
% Step 1: Epoching the data by time window and averaging the epochs
% together at the subject level
mcpa_struct = MCP_to_MCPA(MCP_struct,p.Results.incl_subjects,p.Results.incl_channels,p.Results.time_window,p.Results.baseline_window);

% Step 2: Apply the desired function (e.g., @nanmean) for summarizing time
% window data. You can write custom functions to deal with time- and
% channel-domain data however you want. Default behavior is to apply the
% function along the first dimension of the MCPA pattern, but this can also
% be changed.
mcpa_summ = summarize_MCPA_Struct(p.Results.summary_handle,mcpa_struct);

%% Prep some basic parameters
n_subj = length(p.Results.incl_subjects);
n_sets = size(sets,1);
n_chan = length(p.Results.incl_channels);
try n_cond = length(unique(p.Results.conditions)); catch, n_cond = length(p.Results.conditions); end

%% Set up the results structure which includes a copy of MCPA_pattern

allsubj_results = create_results_struct(mcpa_summ,...
                                        MCP_struct,...
                                        p,...
                                        sets,...
                                        n_subj,...
                                        n_sets,...
                                        n_chan,...
                                        n_cond);

                                                                       
%% Begin the n-fold process: Select one test subj at a time from MCPA struct
for s_idx = 1:length(mcpa_summ.incl_subjects)
    if p.Results.verbose
        fprintf('Running %g feature subsets for Subject %g / %g \n',n_sets,s_idx,n_subj);
    end
    tic;

    % On first fold, initialize the matrix for pairwise results
    if s_idx==1, allsubj_results.accuracy_matrix = nan(n_cond,n_cond,n_sets,n_subj); end
    for set_idx = 1:min(n_sets,p.Results.max_sets)
        
        set_chans = sets(set_idx,:);
        
        % Perform the test for this fold (all possible pairs of conds)
        participant_rsa_matrix = atanh(corr(mcpa_summ.patterns(:,set_chans,s_idx)','rows','pairwise', 'type','Spearman'));
        
        save part_mat participant_rsa_matrix
        
        size(participant_rsa_matrix)
        
        [subj_acc, comparisons] = pairwise_rsa_test(participant_rsa_matrix,semantic_model);

        % Record the results into the results struct
        for comp = 1:size(comparisons,1)
            allsubj_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),set_idx, s_idx) = subj_acc(comp);
        end
    end
        
end
    
    
    %% Progress reporting
    if p.Results.verbose
        fprintf(' %0.1f mins\n',toc/60);
    end


%% Visualization
if p.Results.verbose
    if n_sets > 1 && length(p.Results.conditions)==2
        
        figure
        errorbar(1:size(allsubj_results.accuracy(1).subjXchan,2),mean(allsubj_results.accuracy(1).subjXchan),std(allsubj_results.accuracy(1).subjXchan)/sqrt(size(allsubj_results.accuracy(1).subjXchan,1)),'r')
        hold;
        errorbar(1:size(allsubj_results.accuracy(2).subjXchan,2),mean(allsubj_results.accuracy(2).subjXchan),std(allsubj_results.accuracy(2).subjXchan)/sqrt(size(allsubj_results.accuracy(2).subjXchan,1)),'k')
        title('Decoding Accuracy across all channels: Red = Cond1, Black = Cond2')
        set(gca,'XTick',[1:length(p.Results.incl_channels)])
        set(gca,'XTickLabel',p.Results.incl_channels)
        hold off;
        
        figure
        errorbar(1:size(allsubj_results.accuracy(1).subjXchan,1),mean(allsubj_results.accuracy(1).subjXchan'),repmat(std(mean(allsubj_results.accuracy(1).subjXchan'))/sqrt(size(allsubj_results.accuracy(1).subjXchan,2)),1,size(allsubj_results.accuracy(1).subjXchan,1)),'r')
        hold;
        errorbar(1:size(allsubj_results.accuracy(2).subjXchan,1),mean(allsubj_results.accuracy(2).subjXchan'),repmat(std(mean(allsubj_results.accuracy(2).subjXchan'))/sqrt(size(allsubj_results.accuracy(2).subjXchan,2)),1,size(allsubj_results.accuracy(2).subjXchan,1)),'k')
        title('Decoding Accuracy across all subjects: Red = Cond1, Black = Cond2')
        set(gca,'XTick',[1:p.Results.incl_subjects])
        hold off;
        
    end
end

end
