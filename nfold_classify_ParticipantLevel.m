function allsubj_results = nfold_classify_ParticipantLevel(MCP_struct,varargin)
%% nfold_classify_ParticipantLevel takes an MCP struct and performs
% n-fold cross-validation for n subjects to classify individual
% participants' average response patterns.
%
%
% This wrapper allows the user to choose if
% features will be averaged within-participants to produce a single
% participant-level observation or if individual events will be preserved.
% If the featuers are averaged within-participants, the training set
% is constrained to the number of participants minus 1. Otherwise the
% training set will be (participants - 1) * (number of instances of an
% object).
%
% Several parameters can be changed,
% including which functions are used to generate features and what
% classifier is trained. See Arguments below:
%
% Arguments:
% MCP_struct: either an MCP-formatted struct or the path to a Matlab file
% (.mat or .mcp) containing the MCP_struct.
% incl_features: features to include in the analysis. Default: all features
% incl_subjects: index of participants to include. Default: all participants
% baseline_window: [onset, offset] in seconds. Default [-5,0]
% time_window: [onset, offset] in seconds. Default [2,6]
% conditions: cell array of condition names / trigger #s. Default: {1,2}
% summary_handle: function handle (or char of function name) to specify how
% time-x-feature data should be summarized into features. Default: nanmean
% setsize: number of features to analyze (for subset analyses) Default: all
% test_handle: function handle for classifier. Default: mcpa_classify
% opts_struct: contains additional classifier options. Default: empty struct
% verbose: logical flag to report status updates and results. Default: true
% scale_data: option to perform some for of feature scaling. Defualt: false
% scale_withinSessions: if we do norm data, this allows the user to choose
% if we norm within individual sessions or across a participants session.
% It is only valid to norm across sessions for
% nfold_classify_ParticipantLevel. Default: true
% scale_function: function handle for how to scale the data. Default:
% minMax_scale
% minMax: what min and max to set the data to. Default: [0,1]
% summarize_dimensions: what dimensions and what order to summarize the
% dimensions of mcpa patterns by. Default behavior is to average dimensions
% in the order of this cell array. Default: {'instance', 'time'}


%% Load MCP struct if necessary
if isstring(MCP_struct) || ischar(MCP_struct)
    MCP_struct = load(MCP_struct,'-mat');
    varname = fieldnames(MCP_struct);
    MCP_struct = eval(['MCP_struct.' varname{1}]);
end


%% if data has already been summarized, leave as is. Otherwise, setup MCPA data and summarize it
if ~any(cellfun(@(x) strcmp(x, 'results_struct'), varargin(find(rem(1:length(varargin), 2)))))
    allsubj_results = setup_MCPA_data(MCP_struct,varargin);
    
else
    allsubj_results = varargin{find(rem(1:length(varargin), 2))+1};
end
    
%% Prep some basic parameters
n_subj = length(allsubj_results.incl_subjects);
n_sets = size(allsubj_results.subsets,1);
n_feature = length(allsubj_results.incl_features);
try n_cond = length(unique(allsubj_results.conditions)); catch, n_cond = length(allsubj_results.conditions); end

%% validate classification options
allsubj_results.opts_struct = validate_classification_options_input(MCP_struct, allsubj_results, allsubj_results.suppress_warnings);

%% Folding & Dispatcher
for s_idx = 1:n_subj    
    if allsubj_results.verbose == 1
        fprintf('Running %g feature subsets for Subject %g / %g',n_sets,s_idx,n_subj);
    end
    tic;
    
    %% Split data into test and train 
    % on each fold, one participant's data will be left out as the test
    % set, the rest of the participants data will be combined according to
    % final_dimensions. 

    [train_data, train_labels, test_data, test_labels] = split_test_and_train(s_idx,...
        allsubj_results.conditions,...
        allsubj_results.patterns,...
        allsubj_results.conditions,...
        allsubj_results.final_dimensions,...
        allsubj_results.dimensions, [], []);
    
    % if running this as part of a significance test, permute the group labels
    if allsubj_results.permutation_test
        num_labels = length(train_labels);
        permuted_idx = randperm(num_labels)';
        train_labels = train_labels(permuted_idx);
    end

    
    %% Run feature subsetting search
    for set_idx = 1:min(n_sets,allsubj_results.max_sets)
        %% Progress reporting bit (not important to function. just sanity)
        % Report at every 5% progress
        if allsubj_results.verbose
            status_jump = floor(n_sets/20);
            if ~mod(set_idx,status_jump)
                fprintf(' .')
            end
        end
        % Select the features for this subset
        set_features = allsubj_results.subsets(set_idx,:);
        
        %% Classify
        %inds = pad_dimensions(allsubj_results.final_dimensions, 'feature+time', set_features);
        [predicted_labels, comparisons] = allsubj_results.test_handle(...
                train_data, ...
                train_labels,...
                test_data,...
                test_labels,...
                allsubj_results.opts_struct);
            
        %% Record results 
        if size(predicted_labels,2) > 1 % test labels will be a column vector if we don't do pairwise          
            subj_acc = nanmean(strcmp(predicted_labels(:,1,:), predicted_labels(:,2,:)));
            nan_idx = cellfun(@(x) any(isnan(x)), predicted_labels(:,1,:), 'UniformOutput', false);
            subj_acc(:,:,[nan_idx{1,:,:}]) = nan;

            % Then loop through comparisons and save accuracy to the results struct
            for comp = 1:size(comparisons,1)
                allsubj_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),set_idx, s_idx) = subj_acc(comp);
            end

        else
            subj_acc = double(strcmp(predicted_labels, test_labels));
            nan_idx = cellfun(@isnan, predicted_labels);
            subj_acc(nan_idx) = NaN;
            for cond_idx = 1:n_cond
                cond_acc = nanmean(subj_acc(strcmp(comparisons, test_labels(cond_idx))));
                allsubj_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = cond_acc;
                allsubj_results.accuracy(cond_idx).subjXfeature(s_idx,:) = cond_acc;
            end
        end

       
    end % end set_idx loop
    
    %% Progress reporting
    if allsubj_results.verbose == 1
        fprintf(' %0.1f mins\n',toc/60);
    end
end % end subject loop

%% Visualization
if allsubj_results.verbose
    if n_sets > 1 && length(input_struct.conditions)==2
        
        figure
        errorbar(1:size(allsubj_results.accuracy(1).subjXfeature,2),mean(allsubj_results.accuracy(1).subjXfeature),std(allsubj_results.accuracy(1).subjXfeature)/sqrt(size(allsubj_results.accuracy(1).subjXfeature,1)),'r')
        hold;
        errorbar(1:size(allsubj_results.accuracy(2).subjXfeature,2),mean(allsubj_results.accuracy(2).subjXfeature),std(allsubj_results.accuracy(2).subjXfeature)/sqrt(size(allsubj_results.accuracy(2).subjXfeature,1)),'k')
        title('Decoding Accuracy across all features: Red = Cond1, Black = Cond2')
        set(gca,'XTick',[1:length(input_struct.incl_features)])
        set(gca,'XTickLabel',input_struct.incl_features)
        hold off;
        
        figure
        errorbar(1:size(allsubj_results.accuracy(1).subjXfeature,1),mean(allsubj_results.accuracy(1).subjXfeature'),repmat(std(mean(allsubj_results.accuracy(1).subjXfeature'))/sqrt(size(allsubj_results.accuracy(1).subjXfeature,2)),1,size(allsubj_results.accuracy(1).subjXfeature,1)),'r')
        hold;
        errorbar(1:size(allsubj_results.accuracy(2).subjXfeature,1),mean(allsubj_results.accuracy(2).subjXfeature'),repmat(std(mean(allsubj_results.accuracy(2).subjXfeature'))/sqrt(size(allsubj_results.accuracy(2).subjXfeature,2)),1,size(allsubj_results.accuracy(2).subjXfeature,1)),'k')
        title('Decoding Accuracy across all subjects: Red = Cond1, Black = Cond2')
        set(gca,'XTick',[1:input_struct.incl_subjects])
        hold off;
        
    end
end

end
