function allsubj_results = nfold_classify_ParticipantLevel_permutationTest(results_struct)
%% nfold_classify_ParticipantLevel_permutationTest takes a results struct and performs
% n-fold cross-validation for n subjects to classify individual
% participants' average response patterns, but with permuted labels. 
% All parameters are taken from the results struct


% Arguments:
% results_struct: results struct from previous decoding


%% Prep some basic parameters
n_subj = length(results_struct.incl_subj);
n_sets = size(results_struct.subsets,1);
n_chan = length(results_struct.incl_channels);
try n_cond = length(unique(results_struct.conditions)); catch, n_cond = length(results_struct.conditions); end


%% Begin the n-fold process: Select one test subj at a time from MCPA struct
for s_idx = 1:length(results_struct.incl_subj)
    if results_struct.verbose
        fprintf('Running %g feature subsets for Subject %g / %g',n_sets,s_idx,n_subj);
    end
    tic;
    
    %% Extract training and testing data
    group_subvec = 1:length(results_struct.incl_subj);
    group_subvec(s_idx) = [];
    
    % Set logical flags for indexing the conditions that will be compared.
    % Loop through the whole list of conditions and create flags for each.
    cond_flags = cell(n_cond,1); % These are, for the moment, empty
    
    %% Run over channel subsets
    temp_set_results_cond = nan(n_cond,n_sets,n_chan);
    
    
    
    %% Folding & Dispatcher
    
    if n_cond==2
        [group_data, group_labels, subj_data, subj_labels] = split_data(s_idx,...
                                                                        cond_flags,...
                                                                        results_struct,...
                                                                        n_cond,...
                                                                        results_struct.summed_mcpa.patterns,...
                                                                        group_subvec,...
                                                                        [],...
                                                                        [],...
                                                                        [],...
                                                                        [],...
                                                                        results_struct.summed_mcpa.event_types);
        %% permute the group labels
        num_labels = length(group_labels);
        permuted_idx = randperm(num_labels)';
        group_labels = group_labels(permuted_idx);
        
        %% Run classifier and compare output with correct labels
        for set_idx = 1:min(n_sets,results_struct.max_sets)
            %% Progress reporting bit (not important to function. just sanity)
            % Report at every 5% progress
            if results_struct.verbose
                status_jump = floor(n_sets/20);
                if ~mod(set_idx,status_jump)
                    fprintf(' .')
                end
            end
            
            % Select the channels for this subset
            set_chans = results_struct.subsets(set_idx,:);
           
            
            %% classify
            results_struct.opts_struct = [];
            temp_test_labels = results_struct.test_handle(...
                group_data(:,set_chans), ...
                group_labels,...
                subj_data(:,set_chans),...
                results_struct.opts_struct);
            
            % Compare the labels output by the classifier to the known labels
            temp_acc1 = cellfun(@strcmp,...
                subj_labels(strcmp(strjoin(string(results_struct.conditions{1}),'+'),subj_labels)),... % known labels
                temp_test_labels(strcmp(strjoin(string(results_struct.conditions{1}),'+'),subj_labels))...% classifier labels
                );
            temp_acc2 = cellfun(@strcmp,...
                subj_labels(strcmp(strjoin(string(results_struct.conditions{2}),'+'),subj_labels)),... % known labels
                temp_test_labels(strcmp(strjoin(string(results_struct.conditions{2}),'+'),subj_labels))... % classifier labels
                );
            
            temp_set_results_cond(1,set_idx,set_chans) = nanmean(temp_acc1);
            temp_set_results_cond(2,set_idx,set_chans) = nanmean(temp_acc2);
        end
        
        
        for cond_idx = 1:n_cond
            allsubj_results.accuracy(cond_idx).subsetXsubj(:,s_idx) = nanmean(temp_set_results_cond(cond_idx,:,:),3);
            allsubj_results.accuracy(cond_idx).subjXchan(s_idx,:) = nanmean(temp_set_results_cond(cond_idx,:,:),2);
        end
       
        
    %% Multiple conditions
        
    else
        if s_idx==1, allsubj_results.accuracy_matrix = nan(n_cond,n_cond,n_subj); end
        
        
        num_groups = size(results_struct.summed_mcpa.patterns,1);
        permuted_idx = randperm(num_groups)';

        [subj_acc, comparisons] = pairwise_rsa_test(results_struct.summed_mcpa.patterns(:,:,s_idx),nanmean(results_struct.summed_mcpa.patterns(permuted_idx,:,group_subvec),3));
        

        for comp = 1:size(comparisons,1)
            allsubj_results.accuracy_matrix(comparisons(comp,1),comparisons(comp,2),s_idx) = subj_acc(comp);
        end
    end
    
    if results_struct.verbose
        fprintf(' %0.1f mins\n',toc/60);
    end
end

end








