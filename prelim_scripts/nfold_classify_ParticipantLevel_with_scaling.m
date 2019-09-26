function allsubj_results = nfold_classify_ParticipantLevel(MCP_struct,varargin)
%% nfold_classify_ParticipantLevel takes an MCP struct and performs
% n-fold cross-validation for n subjects to classify individual
% participants' average response patterns. This wrapper assumes that
% features will be averaged within-participants to produce a single
% participant-level observation. Thus the training set is constrained to
% the number of participants minus 1. Several parameters can be changed,
% including which functions are used to generate features and what
% classifier is trained. See Arguments below:
%
% Arguments:
% MCP_struct: either an MCP-formatted struct or the path to a Matlab file
% (.mat or .mcp) containing the MCP_struct.
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
p = inputParser;
addParameter(p,'incl_channels',[1:max(arrayfun(@(x) size(x.fNIRS_Data.Hb_data.Oxy,2),MCP_struct))],@isnumeric);
addParameter(p,'incl_subjects',[1:length(MCP_struct)],@isnumeric);
addParameter(p,'time_window',[2,6],@isnumeric);
addParameter(p,'baseline_window',[-5 0],@isnumeric);
addParameter(p,'conditions',unique(cellstr(char(cellfun(@(x) char(x{:}), arrayfun(@(x) unique({x.Experiment.Conditions.Name},'stable'),MCP_struct, 'UniformOutput',false),'UniformOutput',false))),'stable'),@iscell);
addParameter(p,'summary_handle',@nanmean);
addParameter(p,'setsize',max(arrayfun(@(x) size(x.fNIRS_Data.Hb_data.Oxy,2),MCP_struct)),@isnumeric);
addParameter(p,'max_sets',1000000,@isnumeric);
addParameter(p,'test_handle',@mcpa_classify);
addParameter(p,'opts_struct',[],@isstruct);
addParameter(p,'verbose',true,@islogical);
addParameter(p,'norm_method', 'none', @ischar);

parse(p,varargin{:})

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
n_all_sets = nchoosek(length(p.Results.incl_channels),p.Results.setsize);
size_of_sets_inmem = n_all_sets*p.Results.setsize*8+100;
try
    available_mem = memory;
    available_mem = available_mem.MemAvailableAllArrays;
catch
    try
        [~,w] = unix('free | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        available_mem = (stats(3)+stats(end))*1000;
    catch
        available_mem = 2E9;
    end
end
if size_of_sets_inmem > 0.50*available_mem
    try
        warning('Too many feature sets will cause memory problems. Randomly generating a subset.');
        sets = nan(p.Results.max_sets,p.Results.setsize);
        for i = 1:p.Results.max_sets
            sets(i,:) = randsample(p.Results.incl_channels,p.Results.setsize);
        end
    catch
        error('Too many feature sets will cause memory problems. Reduce setsize or length of incl_channels.');
    end
    
else
    if p.Results.verbose
        fprintf('Generating %g possible sets... ',n_all_sets);
    end
    % Start by generating list of all sets
    tic;
    sets = nchoosek(p.Results.incl_channels,p.Results.setsize);
    if p.Results.verbose
        fprintf('Done\n');
        toc
    end
    % Case where there are more possible sets than the max_sets limit
    if p.Results.max_sets < n_all_sets
        % Randomly select rows to use
        sets_to_choose = randsample(1:size(sets,1),p.Results.max_sets);
        % Select the randomly sampled rows from sets.
        sets = sets(sets_to_choose,:);
        if p.Results.verbose
            fprintf('Selected %g sets to test.\n',p.Results.max_sets);
        end
    end
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
% To-do: write a separate constructor script for this because we need
% results structs from all kinds of testing and should be getting the same
% structure out each time.
allsubj_results = [];
allsubj_results.patterns = mcpa_summ.patterns;
allsubj_results.MCP_data = MCP_struct;
allsubj_results.created = datestr(now);
allsubj_results.test_handle = p.Results.test_handle;
allsubj_results.test_type = 'N-fold (Leave one subject out), Classify participant-level averages';
allsubj_results.setsize = p.Results.setsize;
allsubj_results.func_handle = p.Results.summary_handle;
allsubj_results.incl_channels = mcpa_struct.incl_channels;
allsubj_results.conditions = p.Results.conditions;
allsubj_results.subsets = sets;

for cond_id = 1:n_cond
    allsubj_results.accuracy(cond_id).condition = allsubj_results.conditions(cond_id);
    allsubj_results.accuracy(cond_id).subjXchan = nan(n_subj,n_chan);
    allsubj_results.accuracy(cond_id).subsetXsubj = nan(n_sets,n_subj);
end


%% Begin the n-fold process: Select one test subj at a time from MCPA struct
for s_idx = 1:length(mcpa_summ.incl_subjects)
    if p.Results.verbose
        fprintf('Running %g feature subsets for Subject %g / %g',n_sets,s_idx,n_subj);
    end
    tic;
    %% Extract training and testing data
    group_subvec = 1:length(mcpa_summ.incl_subjects);
    group_subvec(s_idx) = [];
    
    % Set logical flags for indexing the conditions that will be compared.
    % Loop through the whole list of conditions and create flags for each.
    cond_flags = cell(n_cond,1); % These are, for the moment, empty
    group_data = [];
    group_labels = [];
    subj_data = [];
    subj_labels = [];
    
    %% Run over channel subsets
    temp_set_results_cond = nan(n_cond,n_sets,n_chan);
    
    
    %% Folding & Dispatcher: Here's the important part
    % Right now, the data have to be treated differently for 2
    % conditions vs. many conditions. In MCPA this is because 2
    % conditions can only be compared in channel space (or, hopefully,
    % MNI space some day). If there are a sufficient number of
    % conditions (6ish or more), we abstract away from channel space
    % using RSA methods. Then classifier is trained/tested on the RSA
    % structures. This works for our previous MCPA studies, but might
    % not be appropriate for other classifiers (like SVM).
    
    %% Two conditions
    % This block should work with most classifiers, provided adequate
    % data are available. We are making the assumption that
    % subject-level averages are the granularity of data that will be
    % both trained and tested.
    if n_cond==2
        
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
        
        %% Run classifier and compare output with correct labels
        for set_idx = 1:min(n_sets,p.Results.max_sets)
            
            %% Progress reporting bit (not important to function. just sanity)
            % Report at every 5% progress
            if p.Results.verbose
                status_jump = floor(n_sets/20);
                if ~mod(set_idx,status_jump)
                    fprintf(' .')
                end
            end
            % Select the channels for this subset
            set_chans = sets(set_idx,:);
            % norm the data
            
            if strcmp(p.Results.norm_method, 'normalize')  
                mins = min(group_data);
                maxes = max(group_data);
                group_data = (group_data - mins) ./ (maxes - mins);
                subj_data = (subj_data-mins) ./ (maxes-mins);
            elseif strcmp(p.Results.norm_method, 'standardize')
                average = mean(group_data);
                sd = std(group_data);
                group_data = (group_data-average) ./ sd;
                subj_data = (subj_data - average) ./ sd;
            else
                group_data = group_data;
                subj_data = subj_data;
            end
                

            % classify
            
            temp_test_labels = p.Results.test_handle(...
                group_data(:,set_chans), ...
                group_labels,...
                subj_data(:,set_chans),...
                p.Results.opts_struct);
            
            % Compare the labels output by the classifier to the known labels
            temp_acc1 = cellfun(@strcmp,...
                subj_labels(strcmp(strjoin(string(p.Results.conditions{1}),'+'),subj_labels)),... % known labels
                temp_test_labels(strcmp(strjoin(string(p.Results.conditions{1}),'+'),subj_labels))...% classifier labels
                );
            temp_acc2 = cellfun(@strcmp,...
                subj_labels(strcmp(strjoin(string(p.Results.conditions{2}),'+'),subj_labels)),... % known labels
                temp_test_labels(strcmp(strjoin(string(p.Results.conditions{2}),'+'),subj_labels))... % classifier labels
                );
            
            % Temporary results from each set are stored in a n_sets x n_chan
            % matrix, so that averaging can be done both across sets (to
            % determine channel mean performance) and across channels (to
            % determine set mean performance)
            temp_set_results_cond(1,set_idx,set_chans) = nanmean(temp_acc1);
            temp_set_results_cond(2,set_idx,set_chans) = nanmean(temp_acc2);
        end
        for cond_idx = 1:n_cond
            allsubj_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
            allsubj_results.accuracy(cond_idx).subjXchan(s_idx,:) = nanmean(temp_set_results_cond(cond_idx,:,:),2);
        end
        
        %% Multiple conditions
        % A bit of complication for how this block should run. If we want
        % an RSA-based classifier, we can either do the all-possible-2way
        % comparisons approach OR we can try doing structural alignment of
        % the whole test dataset (similarity structure) to the training
        % dataset (another similarity structure), a la Zinszer et al.,
        % 2016, Journal of Cognitive Neuroscience (fMRI-based translation).
        %
        % If we don't want to do RSA based (i.e., stay in channel or MNI
        % space), then we need to ask whether we're doing all-possible-2way
        % comparisons or some n-alternative-forced-choice test with chance
        % performance at 1/n.
        %
        % No graceful way to handle these branching decisions yet.  We are
        % also still making the assumption that subject-level averages are
        % the granularity of data that will be both trained and tested.
    else
        % TO DO: Write the multiclass dispatcher here
        %
        % For now, just adapting the Neurophotonics script which has
        % result-writing built into it. This is not a good long term
        % solution because it breaks the modularity of the software
        % (and does nothing to support n-fold for all the other
        % possible classifiers we might want.
        %allsubj_results = pairwise_rsa_leaveoneout(mcpa_summ.patterns);
        
        % On first fold, initialize the matrix for pairwise results
        if s_idx==1, allsubj_results.accuracy_matrix = nan(n_cond,n_cond,n_subj); end
        
        % Perform the test for this fold (all possible pairs of conds)
        [subj_acc, comparisons] = pairwise_rsa_test(mcpa_summ.patterns(:,:,s_idx),nanmean(mcpa_summ.patterns(:,:,group_subvec),3));
        
        % Record the results into the results struct
        for comp = 1:size(comparisons,1)
            allsubj_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),s_idx) = subj_acc(comp);
        end
        for cond_idx = 1:n_cond
            %mean_cond_acc = nanmean([allsubj_results.accuracy_matrix(cond_idx,:,s_idx)';allsubj_results.accuracy_matrix(:,cond_idx,s_idx)]);
            %allsubj_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
            %allsubj_results.accuracy(cond_idx).subjXchan(s_idx,:) = nanmean(temp_set_results_cond(cond_idx,:,:),2);
        end
        
    end
    
    
    %% Progress reporting
    if p.Results.verbose
        fprintf(' %0.1f mins\n',toc/60);
    end
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