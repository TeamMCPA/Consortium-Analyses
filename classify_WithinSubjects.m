function results_struct = classify_WithinSubjects(results_struct)
%% nfold_classify_ParticipantLevel takes in a struct containing the MCPA summarized patterns 
% and classification parameters parsed from cross_validate,
% nest_cross_validate, or significance_test. It will then perform within
% subjects cross validation. It can handle leave-one-session out and k-Fold
% CV.
% It then returns the original struct with the accuracy of each
% fold. 

% created by Anna Herbolzheimer and Ben Zinszer Fall 2020

%% Prep some basic parameters
n_subj = length(results_struct.incl_subjects);
n_sets = size(results_struct.subsets,1);
n_feature = length(results_struct.incl_features);
s = size(results_struct.patterns);
n_sessions = s(end-1);
try n_cond = length(unique(results_struct.conditions)); catch, n_cond = length(results_struct.conditions); end

for s_idx = 1:n_subj
    
    if results_struct.verbose
        fprintf('Running %g feature subsets for Subject %g / %g \n',n_sets,s_idx,n_subj);
    end
    tic;
    
    inds = pad_dimensions(results_struct.dimensions, 'subject', s_idx);
    subject_patterns = results_struct.patterns(inds{:});
    
    %% set up indexing for either k-fold(kf) or leave one out (loo) CV
    %% K-Fold (kf)
    if strcmp(results_struct.approach, 'kf')   
        subject_labels = repmat(results_struct.event_types, (size(subject_patterns,1)/8),1);

        % remove empty rows
        remove = [];
        for r = 1:size(subject_patterns,1)
            if sum(isnan(subject_patterns(r,:))) > 0
               remove = [remove; r];
            end
        end

        subject_patterns(remove,:) = [];
        subject_labels(remove) = [];

        % find indices for each fold
        num_data = size(subject_patterns,1);
        num_in_fold = floor(num_data*results_struct.test_percent);
        num_folds = floor(num_data/num_in_fold);
        fold_end_idx_array = [num_in_fold: num_in_fold: num_in_fold*num_folds];
        fold_end_idx_array(end) = num_data;
        fold_start_idx_array = [1:num_in_fold: num_in_fold*num_folds];

        
        % should we randomize the rows?
        if strcmp(results_struct.randomized_or_notrand, 'randomized')
            rng('default')
            rand_inds = randperm(size(subject_patterns,1));
            subject_patterns(rand_inds, :);
            subject_labels(rand_inds);
        end
        
        
        % should we balance the classes?
        if results_struct.balance_classes
            num_repetitions = num_data/length(results_struct.conditions);
            new_kfold_mat = validate_balanced_classes(fold_start_idx_array,...
                fold_end_idx_array,...
                results_struct,...
                num_folds,...
                subject_labels,...
                num_repetitions);
        end
        
        % define where test and train labels are derived from
        event_types = subject_labels;

    %% Leave-One-Out (loo)
    elseif strcmp(results_struct.approach, 'loo')   
        % randomize trials, if needed
        if strcmp(results_struct.randomized_or_notrand, 'randomized') 
            new_index = randperm(size(subject_patterns,ndims(subject_patterns)));
            subject_patterns = subject_patterns(:,:,new_index);
        end
        
        % define where test and train labels are derived from
        event_types = results_struct.event_types;
        
        % define the fold
        num_folds = size(subject_patterns, ndims(subject_patterns));
    end
    
    %% Begin cross-validation
    for folding_idx = 1:num_folds         
        %% Define fold_idx: indices in subject_patterns that will be test data
        temp_set_results_cond = nan(n_cond,n_sets,n_feature);
        
        if strcmp(results_struct.approach, 'kf')
            fold = new_kfold_mat(folding_idx,~isnan(new_kfold_mat(folding_idx,:)));
        else
            fold = folding_idx;
        end

        %% Split data into test and train 
        % on each fold, one participant's data will be left out as the test
        % set, the rest of the participants data will be combined according to
        % final_dimensions. 
        
        [train_data, train_labels, test_data, test_labels] = split_test_and_train(fold,...
            results_struct.conditions,...
            subject_patterns,... 
            event_types,...
            {},...
            results_struct.dimensions, [], []);
                
        % Skip an iteration if all of train_data or all of test_data is NaN
        if sum(~isnan(train_data(:)))==0 || sum(~isnan(test_data(:)))==0
            fprintf('skip sub %d fold %d \n', s_idx, folding_idx);
            continue; 
        end
                   
        %% Take a random subset, if needed
        
        % Only for K-fold. If what you're trying to sample is greater
        % than the available dataset, this function will sample with
        % replacement within training data and within testing data separately,
        % so that no single data point is used in both training and testing
        % datasets. If what you'er trying to sample is less than the
        % available dataset, this function will sample without replacement.
        
        if strcmp(results_struct.approach, 'kf') && (~isfield(results_struct, 'randomsubset')||isempty(results_struct.randomsubset))
            allsubj_results.kf_randomsubset = 'no: using the whole dataset';
            allsubj_results.kf_sampling = 'no sampling';  
            
        elseif strcmp(results_struct.approach, 'kf') && isfield(results_struct, 'randomsubset')
            n_test = size(test_data,1);
            n_train = size(train_data,1);
            
            subset_test = floor(results_struct.randomsubset*results_struct.test_percent);
            subset_train = floor(results_struct.randomsubset*(1-results_struct.test_percent));
            allsubj_results.kf_randomsubset = sprintf('yes: %d out of %d test data and %d out of %d train data', subset_test, n_test, subset_train, n_train); 

            if subset_test <= n_test
                sampled_test_data = randsample(1:n_test,subset_test,false);
                sampled_train_data = randsample(1:n_train,subset_train,false);
                
                test_data = test_data(sampled_test_data,:);
                test_labels = test_labels(sampled_test_data);
                train_data = train_data(sampled_train_data,:);
                train_labels = train_labels(sampled_train_data);
                allsubj_results.kf_sampling = 'sample without replacement';  
                
            elseif  subset_test > n_test 
                sampled_test_data = randsample(1:n_test,subset_test,true);
                sampled_train_data = randsample(1:n_train,subset_train,true);
                
                test_data = test_data(sampled_test_data,:);
                test_labels = test_labels(sampled_test_data);
                train_data = train_data(sampled_train_data,:);
                train_labels = train_labels(sampled_train_data);
                
                allsubj_results.kf_sampling = 'sample with replacement';   
            end
            
        end

        %% Run classifier and compare output with correct labels        
        for set_idx = 1:min(n_sets,results_struct.max_sets)    

            if results_struct.verbose
                status_jump = floor(n_sets/20);
                if ~mod(set_idx,status_jump)
                    fprintf(' .')
                end
            end
            
            % Select the features for this subset
            set_features = results_struct.subsets(set_idx,:);

            %% Classify
            inds = pad_dimensions(results_struct.final_dimensions, 'feature', set_features);
            [predicted_labels, comparisons] = results_struct.test_handle(...
                    train_data(inds{:}), ...
                    train_labels,...
                    test_data(inds{:}),...
                    test_labels,...
                    results_struct.opts_struct);

            %% Record results 
            results_struct = save_accuracy(results_struct,...
                comparisons,...
                predicted_labels,...
                set_idx,...
                fold,... 
                s_idx);

        end % end set_idx
        
    if results_struct.verbose
        fprintf(' %0.1f mins\n',toc/60);
    end
    end % end session loop
    
end % end subject loop


end