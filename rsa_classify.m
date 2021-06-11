function [classification, comparisons] = rsa_classify(model_data, model_labels, test_data, test_labels, opts)
%% rsa_classify implements a correlation-based, similarity-space classifier
% following Zinszer, Bayet, Emberson, Raizada, & Aslin's (2018, Neurophotonics)
% method for all-possible-pairwise-comparisons and Zinszer, Anderson, Kang,
% Wheatley, & Raizada's (2016, JoCN) approach for n-way comparison
%
% opts: a struct that contains options for the classifier.
%
% opts.comparison_type determines how the features will be abstracted:
%    'comparison_type' if true, will use correlation-based similarity
%    space (larger values more similar). opts.metric can be set to 
%    'spearman', 'pearson', or 'kendall'.
%    'comparison_type' if false, will use pdist-based dissimilarity matrix
%    (larger values are less similar). Any of the distance metrics from the
%    pdist function (e.g., 'cosine', 'euclidean', 'seuclidean') can be set
%    for the opts.metric value.
%    Default setting is opts.similiarty_space=true, opts.metric='spearman'
% 
%    If similarity space is used, values are adjusted with Fisher's r-to-z
%    transform (the hyperbolic arctangent of the correlation coefficient).
%
% Only opts.exclusive=true mode is currently available.
%
% opts.pairwise if true will make all possible pairwise contrasts between
% the classes. If false, it will test all n-way comparison of class labels
% for permutation of labels with highest correlation to the training set.
%
% If opts.tiebreak is set to true (default), the correlation coefficients
% are randomly adjusted by <1% of the smallest observed difference to
% prevent exact matches and thus prevent ties in the classification. If
% opts.tiebreak is set to false, the classifier will prefer classes
% appearing earlier in the training set.
%
% If opts.verbose is set to true (default is false), Fisher-adjusted
% (hyperbolic arctangent) correlation matrices are displayed for all
% subjects, sessions, etc.


%% Pull a list of all the unique classes / conditions, preserving order
model_classes = unique(model_labels(:),'stable');

%% check to see the orders of test and train data - this is to see if they match
% get number of categories

[~,train_order] = sort(model_labels);
[~,test_order] = sort(test_labels);

% sort test data based on the re-ordered data
test_labs = test_labels;%(test_order);
test_dat = test_data;%(test_order, :,:,:);

model_labs = model_labels;%(train_order);
model_dat = model_data;%(train_order, :,:,:);

%% remove where we have all NaNs
% 
% empty_x_vals_model = [];
% empty_y_vals_model = [];
% for i = 1:size(model_data,4)
%     for j = 1:size(model_data,3)
%         [x,y] = find(isnan(model_dat(:,:,j,i)));
%         if length(unique(x)) == size(model_dat,1) && length(unique(y)) == size(model_dat,2)
%             continue;
%         else
%             empty_x_vals_model = [empty_x_vals_model; x];
%             empty_y_vals_model = [empty_y_vals_model; y];
%         end
%     end      
% end
% empty_x_vals_model = unique(empty_x_vals_model);
% empty_y_vals_model = unique(empty_y_vals_model);
% 
% 
% empty_x_vals_test = [];
% empty_y_vals_test = [];
% for i = 1:size(test_data,4)
%     for j = 1:size(test_data,3)
%         [x,y] = find(isnan(test_data(:,:,j,i)));
%         
%         
%         if length(unique(x)) == size(test_data,1) && length(unique(y)) == size(test_data,2)
%             continue;
%         else
%             empty_x_vals_test = [empty_x_vals_test; x];
%             empty_y_vals_test = [empty_y_vals_test; y];
%         end
%     end      
% end
% 
% empty_x_vals_test = unique(empty_x_vals_test);
% empty_y_vals_test = unique(empty_y_vals_test);
% 
% cols = 1:size(model_data,2);
% rows = 1:size(model_data,1);
% 
% if size(model_data,2) > 8
%     if length(empty_x_vals_model) == size(model_data,1) && length(empty_y_vals_model) ~= size(model_data,2)
%         remove columns
%         remove_cols = union(empty_y_vals_model, empty_y_vals_test);
%         keep = ~ismember(cols, remove_cols);
%         model_dat = model_dat(:,keep', :,:);
%         test_dat = test_dat(:, keep', :,:);   
%     elseif length(empty_x_vals_model) == size(model_data,1) && length(empty_y_vals_model) == size(model_data,2)
%         this is an empty session
%         warning('There is no data in this session')
%     end
% else
%     if length(empty_x_vals_test) == size(model_data,1) && length(empty_y_vals_test) ~= size(test_data,2)
%         test_cols = 1:size(test_dat,2);
%         remove_cols = union(empty_y_vals_model, empty_y_vals_test);
%         keep = ~ismember(test_cols, remove_cols);
%         test_dat = test_dat(:, keep', :,:);
%     end 
% end
    
%% Build similarity structures
% Transform the model_data into similiarty or dissimilarity structures for each session by
% correlating between conditions or finding pairwise differences between conditions and then 
% average the similarity/dissimilarity structures together into a single model.

if strcmp(opts.comparison_type, 'correlation') % create similarity structures
    %  iterate through all the layers (3rd dimension) and create
    % correlation matrices 
    
    if ndims(model_dat) ~= 2 && size(model_dat,1) == length(model_classes)
        % if we're doing leave one out (participant or session), first want
        % to correlate the model data within session, then average across sessions
        model_mat = nan(size(model_dat,1),size(model_dat,1),size(model_dat,3),size(model_dat,4));
        for i = 1: (size(model_dat,3)*size(model_dat,4))
            model_mat(:,:,i) = corr(model_dat(:,:,i)','rows','pairwise', 'type', opts.metric);
        end
        
        
        % then average across each participants' sessions
        training_matrix = nanmean(model_mat,3);
        % then average across participants
        training_matrix = nanmean(training_matrix,4);
        training_matrix = atanh(training_matrix);
        
        
        % then repeat to create test data
        test_mat = nan(size(test_dat,1),size(test_dat,1),size(test_dat,3),size(test_dat,4));
        for i = 1: (size(test_dat,3)*size(test_dat,4))
            test_mat(:,:,i) = corr(test_dat(:,:,i)','rows','pairwise', 'type', opts.metric);
        end
        test_matrix = nanmean(test_mat,3);
        test_matrix = nanmean(test_matrix,4);

        test_matrix = atanh(test_matrix);
        
    elseif ndims(model_dat) == 2 && size(model_dat,1) ~= length(model_classes)
        % if doing kfold for withinsubjects cross validation, first need to aggregate over repetitions for each category  
        temp_model_dat = nan(length(model_classes), size(model_dat,2));
        temp_test_dat = nan(length(model_classes), size(test_dat,2));
        for cl = 1:length(model_classes)
            model_dat_for_this_class = model_dat(strcmp(model_labs, model_classes{cl}),:);
            temp_model_dat(cl,:) = nanmean(model_dat_for_this_class,1);

            test_dat_for_this_class = test_dat(strcmp(test_labs, model_classes{cl}),:);
            temp_test_dat(cl,:) = nanmean(test_dat_for_this_class,1);
        end

        % then create test and train matrices in similarity space
        test_matrix = atanh(corr(temp_test_dat'));
        training_matrix = atanh(corr(temp_model_dat'));

        model_labs = model_classes;
        test_labs = model_classes;
    else
        % if using a premade training model (like a semantic model), don't
        % need to do anything to model_dat
        training_matrix = model_dat;
        
        % then create test data
        test_mat = nan(size(test_dat,1),size(test_dat,1),size(test_dat,3),size(test_dat,4));
        for i = 1: (size(test_dat,3)*size(test_dat,4))
            test_mat(:,:,i) = corr(test_dat(:,:,i)', 'type', opts.metric);
        end
        test_matrix = nanmean(test_mat,3);
        test_matrix = nanmean(test_matrix,4);

        test_matrix = atanh(test_matrix);
        
    end
else % create dissimilarity structures
    
    if ndims(model_data) > 2
         % if we're doing leave one out (participant or session), first want
        % to loop through each participants' sessions to find pairwise differences
        model_mat = nan(size(model_dat,1),size(model_dat,1),size(model_dat,3),size(model_dat,4));
        for i = 1: (size(model_dat,3)*size(model_dat,4))
            model_mat(:,:,i) = squareform(pdist(model_dat(:,:,i), opts.metric));
        end
        % then average across each participants' sessions
        training_matrix = nanmean(model_mat,3);
        % then average across participants
        training_matrix = nanmean(training_matrix,4);
        
        % then repeat to create test data
        test_mat = nan(size(test_dat,1),size(test_dat,1),size(test_dat,3),size(test_dat,4));
        for i = 1: (size(test_dat,3)*size(test_dat,4))
            test_mat(:,:,i) = squareform(pdist(test_dat(:,:,i), opts.metric));
        end
        test_matrix = nanmean(test_mat,3);
        test_matrix = nanmean(test_matrix,4);
        
    elseif ndims(model_dat) ~= 2 && size(model_dat,1) ~= length(model_classes)
        % if doing kfold for withinsubjects cross validation, first need to aggregate over repetitions for each category  
        temp_model_dat = nan(length(model_classes), size(model_dat,2));
        temp_test_dat = nan(length(model_classes), size(test_dat,2));
        for cl = 1:length(model_classes)
            model_dat_for_this_class = model_dat(strcmp(model_labs, model_classes{cl}),:);
            temp_model_dat(cl,:) = nanmean(model_dat_for_this_class,1);

            test_dat_for_this_class = test_dat(strcmp(test_labs, model_classes{cl}),:);
            temp_test_dat(cl,:) = nanmean(test_dat_for_this_class,1);
        end

        % then create dissimilarity matrices
        test_matrix = squareform(pdist(temp_test_dat, opts.metric));
        training_matrix = squareform(pdist(temp_model_dat, opts.metric));

        model_labs = model_classes;
        test_labs = model_classes;
    else
        % if using a premade training model (like a semantic model), don't
        % need to do anything to model_dat
        training_matrix = model_dat;
        
        % then repeat to create test data
        test_mat = nan(size(test_dat,1),size(test_dat,1),size(test_dat,3),size(test_dat,4));
        for i = 1: (size(test_dat,3)*size(test_dat,4))
            test_mat(:,:,i) = squareform(pdist(test_dat(:,:,i), opts.metric));
        end
        test_matrix = nanmean(test_mat,3);
        test_matrix = nanmean(test_matrix,4);
    end
    
   
end


%% Visualize the matrices

if opts.verbose > 1
    
    plot_idx = 1;
    figure
    for session_idx = 1:size(model_data,3)
        for subject_idx = 1:size(model_data,4)
            subplot(size(model_data,3),size(model_data,4),plot_idx);
            imagesc(model_mat(:,:,session_idx, subject_idx))
            title(['Subj ' num2str(subject_idx) ' Sess ' num2str(session_idx)])
            xticklabels([])
            yticklabels([])
            caxis([-.5,.5])
            colorbar('hot')
            plot_idx = plot_idx + 1;
        end
    end
    
    plot_idx = 1;
    figure
    for session_idx = 1:size(test_data,3)
        for subject_idx = 1:size(test_data,4)
            subplot(size(test_data,3),size(test_data,4),plot_idx);
            imagesc(test_mat(:,:,session_idx, subject_idx))
            title(['Subject ' num2str(subject_idx) ' Session ' num2str(session_idx)])
            xticklabels([])
            yticklabels([])
            caxis([-.5,.5])
            colorbar('hot')
            plot_idx = plot_idx + 1;
        end
    end
    
end
if opts.verbose
    figure();
    % Training data
    subplot(1,2,1)
    imagesc(training_matrix)
    title('Training Data')
    xticklabels(model_labels)
    yticklabels(model_labels)
    %caxis([min(min(tril(training_matrix,-1))),max(max(tril(training_matrix,-1)))])
    caxis([-.5,.5])
    
    colorbar('hot')
    [i, j, ~] = find(~isnan(training_matrix));
    text(i-.3,j,num2str(round(training_matrix(~isnan(training_matrix)),2)));
    
    % Test data
    subplot(1,2,2)
    imagesc(test_matrix)
    title('Test Data')
    xticklabels(test_labels)
    yticklabels(test_labels)
    %caxis([min(min(tril(test_matrix,-1))),max(max(tril(test_matrix,-1)))])
    caxis([-.5,.5])
    colorbar('hot')
    [i, j, ~] = find(~isnan(test_matrix));
    text(i-.3,j,num2str(round(test_matrix(~isnan(test_matrix)),2)));
end

%% Sanity Check
if sum(isnan(test_matrix(:)))==(numel(test_matrix)-8) || sum(isnan(test_matrix(:)))==numel(test_matrix) || sum(isnan(training_matrix(:)))==numel(training_matrix)
     warning('One or both input matrices contains all NaN values. this one triggered');
end

%% Save out the classification results based on greatest correlation coefficient for each test pattern
% Initialize empty cell matrix for classifications

if ~isfield(opts,'pairwise') || ~opts.pairwise
    
    % Number of stimulus classes to compare, based on the number of rows
    number_classes = size(test_matrix,1);
    
    % Too many permutations if more than 10 classes. Throw error and leave.
    if number_classes>10, error('n-way classification does not work for more than 10 classes. Use opts.pairwise=true'); end
    
    % Generate a list of every pairwise comparison and the results array
    list_of_comparisons = perms(1:number_classes);
    number_of_comparisons = size(list_of_comparisons,1);
    results_of_comparisons = nan(number_of_comparisons,1);
    
    % Convert training set to a vector for the comparisons
    train_vec = training_matrix(logical(tril(ones(size(training_matrix)),-1)));
    
    % Iterate through the different Permutations of test
    for perm_idx = 1:number_of_comparisons
        tmp_test_matrix = test_matrix(list_of_comparisons(perm_idx,:),list_of_comparisons(perm_idx,:));
        tmp_test_vec = tmp_test_matrix(logical(tril(ones(size(tmp_test_matrix)),-1)));
        results_of_comparisons(perm_idx) = corr(train_vec,tmp_test_vec,'rows','pairwise');
    end
    
    % Choose the best of all the permutations of labels
    [rating, best_perm] = max(results_of_comparisons);
    classification = model_classes(list_of_comparisons(best_perm,:));
    
    %accuracy = strcmp(classification,test_labels);
    comparisons = test_labels';
    
    % put the labels back in the order they were put in as
    [~,reorder_test] = sort(test_order);
    classification = classification(reorder_test);
    comparisons = reorder_test;
    
    
else
    
    [accuracy, comparisons] = pairwise_rsa_test(test_matrix,training_matrix);
    
    classification = comparisons;
    
    for acc = 1:length(accuracy)
        if isnan(accuracy(acc))
            classification(acc,:) = nan;
        else
            if ~accuracy(acc)
                classification(acc,:) = classification(acc,end:-1:1);
            end
        end
    end
   
    classification = num2cell(classification);
    for class = 1:size(classification,1)
        if ~any(isnan([classification{class,:}]))
            classification(class,:) = model_classes([classification{class,:}]);
        end
    end
    
    true_labels = model_classes(comparisons); 
    results_of_comparisons = cell(size(classification,2), 2, size(classification,1));
    for comp = 1:size(classification,1)
        results_of_comparisons(:,1,comp) = classification(comp,:)';
        results_of_comparisons(:,2,comp) = true_labels(comp,:)';
    end
    classification = results_of_comparisons;
    
    
end


end
