function allsubj_results = nfold_classify_WithinSubjects(MCP_struct, varargin)
%% function to classify sessions from individual participants

%% Parse out the input data
p = parse_inputs(MCP_struct, varargin{:});

%% Setting up the combinations of channel subsets
unmapped_sets = find_sets(p.Results);
sets = map_values(p, unmapped_sets);

%% norm check - do we want to scale individual participant data?
if p.Results.norm_data
    MCP_struct = scale_individuals(MCP_struct, p.Results);
end

%% Build MCPA struct for all subjects in the MCP
mcpa_struct = MCP_to_MCPA(MCP_struct,...
                         p.Results.incl_subjects,...
                         p.Results.incl_channels,...
                         p.Results.time_window,...
                         p.Results.baseline_window);

%% decide how we want to concatenate or average over our dimensions
% first see if the user specified the summarizing dimensions. If not,
% recommend what dimensions to average over
p.Results
if ~isfield(p.Results, 'summarize_dimensions') || ~isempty(p.Results.summarize_dimensions)
    summarize_dimensions = p.Results.summarize_dimensions;
else
    isWithinSubjects = true;
    summarize_dimensions = recommend_summarize_dimensions(p.Results, isWithinSubjects);
end

%% summarize the MCPA struct
mcpa_summ = summarize_MCPA_Struct(p.Results.summary_handle,...
                                    mcpa_struct,...
                                    summarize_dimensions);

%% Prep some basic parameters
n_subj = length(p.Results.incl_subjects);
n_sets = size(sets,1);
n_chan = length(p.Results.incl_channels);
s = size(mcpa_summ.patterns);
n_sessions = s(end-1);
try n_cond = length(unique(p.Results.conditions)); catch, n_cond = length(p.Results.conditions); end

%% Set up the results structure which includes a copy of MCPA_pattern
allsubj_results = create_results_struct(true,...
                                      mcpa_summ,...
                                      MCP_struct,...
                                      p,...
                                      sets,...
                                      n_subj,...
                                      n_sets,...
                                      n_chan,...
                                      n_cond);
                                          

func = dbstack;
current_folding_function = func.name;
allsubj_results.test_type = current_folding_function;

%% now begin the fold

for s_idx = 1:length(MCP_struct)
    
    if p.Results.verbose
        fprintf('Running %g feature subsets for Subject %g / %g \n',n_sets,s_idx,n_subj);
    end
    tic;
    
    %%
    if length(size(mcpa_summ.patterns)) == 5
        subject_pattern = mcpa_summ.patterns(:,:,:,:,s_idx);
    else
        subject_pattern = mcpa_summ.patterns(:,:,:,s_idx);
    end
    
    for session_idx = 1:length(MCP_struct(s_idx).Experiment.Runs)
        %% Extract training and testing data
        group_subvec = 1:length(MCP_struct(s_idx).Experiment.Runs);
        group_subvec(session_idx) = [];

        % Set logical flags for indexing the conditions that will be compared.
        % Loop through the whole list of conditions and create flags for each.
        cond_flags = cell(n_cond,1); % These are, for the moment, empty
        
        %% Run over channel subsets
        
        if n_cond == 2
            group_data = [];
            group_labels = [];
            subj_data = [];
            subj_labels = [];
            temp_set_results_cond = nan(n_cond,n_sets,n_chan);

            
            [group_data, group_labels, subj_data, subj_labels] = split_data(session_idx,...
                                                                                cond_flags,...
                                                                                p.Results,...
                                                                                n_cond,...
                                                                                subject_pattern,...
                                                                                group_subvec,...
                                                                                group_data,...
                                                                                group_labels,...
                                                                                subj_data,...
                                                                                subj_labels,...
                                                                                mcpa_summ.event_types);
                                                                            

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

                %% classify

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
                allsubj_results.accuracy(cond_idx).subjXsession(s_idx,session_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
            end

        else
            
            [group_data, ~, subj_data, ~] = split_data_clean(session_idx,...
                                                            cond_flags,...
                                                            p.Results,...
                                                            n_cond,...
                                                            subject_patterns,...
                                                            group_subvec,...
                                                            [],...
                                                            [],...
                                                            [],...
                                                            [],...
                                                            mcpa_summ.event_types);
                                                        
                                                        
                                                        
            
            if strcmp(func2str(p.Results.test_handle),'pairwise_rsa_test') 
                num_sets = min(n_sets,p.Results.max_sets);
                
                if s_idx == 1 && session_idx == 1, allsubj_results.accuracy_matrix = nan(n_cond,n_cond,num_sets, n_sessions, n_subj); end
                
                for set_idx = 1:num_sets
                    % now do classification for each channel subset
                    set_chans = sets(set_idx,:);
                    
                    [subj_acc, comparisons] = pairwise_rsa_test(subj_data(:,set_chans),...
                                                                group_data(:,set_chans));
                  
                    for comp = 1:size(comparisons,1)
                        
                        allsubj_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),set_idx,session_idx,s_idx) = subj_acc(comp);
                    end
                end
            else
                [group_data, group_labels, subj_data, subj_labels] = split_data_clean(session_idx,...
                                                            cond_flags,...
                                                            p.Results,...
                                                            n_cond,...
                                                            subject_patterns,...
                                                            group_subvec,...
                                                            [],...
                                                            [],...
                                                            [],...
                                                            [],...
                                                            mcpa_summ.event_types);
                                                        
                for set_idx = 1:min(n_sets,p.Results.max_sets)
                    
                    set_chans = sets(set_idx,:);
                    
                    temp_test_labels = p.Results.test_handle(...
                    group_data(:,set_chans), ...
                    group_labels,...
                    subj_data(:,set_chans),...
                    p.Results.opts_struct);

                    % Temporary results from each set are stored in a n_sets x n_chan
                    % matrix, so that averaging can be done both across sets (to
                    % determine channel mean performance) and across channels (to
                    % determine set mean performance)

                    
                    for cond_idx = 1:n_cond
                        temp_acc = cellfun(@strcmp,...
                        subj_labels(strcmp(strjoin(string(p.Results.conditions{cond_idx}),'+'),subj_labels)),... % known labels
                        temp_test_labels(strcmp(strjoin(string(p.Results.conditions{cond_idx}),'+'),subj_labels))...% classifier labels
                        );
                    
                        temp_set_results_cond(cond_idx,set_idx,set_chans) = nanmean(temp_acc);
                    
                        allsubj_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
                        allsubj_results.accuracy(cond_idx).subjXchan(s_idx,:) = nanmean(temp_set_results_cond(cond_idx,:,:),2);
                        allsubj_results.accuracy(cond_idx).subjXsession(s_idx,session_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
                    end
                end
            end            
        end
    end
    
end


end
            

            
