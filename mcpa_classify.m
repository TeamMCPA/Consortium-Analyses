function classification = mcpa_classify(model_data, model_labels, test_data, opts)
%% mcpa_classify implements a correlation-based, feature-space classifier
% following Emberson, Zinszer, Raizada & Aslin's (2017, PLoS One) method.
%
% Currently implementation only works for two conditions at a time.
% Multiple conditions is under developments (fingers crossed...)

%% If the options struct is not provided, set default parameters
if ~exist('opts','var') || isempty(opts)
    opts = struct;
    opts.corr_stat = 'spearman';
    opts.exclusive = false;
end

%% Average across training data to get model features for each class
model_classes = unique(model_labels,'stable');
model_pattern1 = nanmean(model_data(strcmp(model_classes{1},model_labels),:),1)';
model_pattern2 = nanmean(model_data(strcmp(model_classes{2},model_labels),:),1)';

%% Perform correlation (default: Pearson) between two model patterns and the test patterns
% This is a quick and easy way to compute all the test items at once.

% 1. Build a matrix with [ModelA, ModelB, Test1, Test2,... TestN]
% 2. Run correlation over all of them, and get a the matrix of corr coeffs.
% 3. First column is Model A vs. all others, Second column is Model B vs.
%	all others.
corr_matrix = atanh(corr([model_pattern1,model_pattern2,test_data'],'type',opts.corr_stat,'rows','pairwise'));

% Isolate the first and second columns (Matrix A and Matrix B), and the
% rows 3 to end (correlations to Matrix A and B of each test vector). Thus
% each row of column 1 is Matrix A vs. Test r (for row r). Each row of
% column 2 is Matrix B vs. Test r (for row r).
test_model_corrs = corr_matrix(3:end,1:2);

%% Save out the classification results based on greatest correlation coefficient for each test pattern
% Initialize empty cell matrix for classifications
classification = cell(size(test_data,1),1);

if opts.exclusive && length(model_classes)==size(test_data,1)
    
    %% Fill in the case for exclusive labels
    % wherein each class can only be assigned to one row of test data
    % (e.g., at Participant level when you are looking at
    % condition-averaged data, and only one pattern per condition)
    if size(test_model_corrs,1)==2
        if trace(test_model_corrs) > trace(rot90(test_model_corrs))
            classification(1) = model_classes(1);
            classification(2) = model_classes(2);
        elseif trace(test_model_corrs) < trace(rot90(test_model_corrs))
            classification(1) = model_classes(2);
            classification(2) = model_classes(1);
        else
            % If both options are equal, randomly assign the two labels to 
            % the two observations.
            order = randperm(2); % returns [1 2] or [2 1] with equal probability
            classification(1) = model_classes(order(1));
            classification(2) = model_classes(order(2));
        end
    else
        % The search-all-label-permutations method would work here for 3 to
        % 10 classes, but the search space becomes too large after 10.
        disp('Currently no method for exclusive labeling with >2 test cases');
    end
    
else
    
    %% Otherwise, just label by best fit for each test_data row
    % Classify based on greatest correlation coefficient
    classification(test_model_corrs(:,1)>test_model_corrs(:,2)) = model_classes(1);
    classification(test_model_corrs(:,1)<test_model_corrs(:,2)) = model_classes(2);
    
    % If correlations to each model pattern are equal, enter NaN
    % temporarily and then replace with random labels
    classification(test_model_corrs(:,1)==test_model_corrs(:,2)) = {NaN};
    num_nans = sum(isnan(classification));
    if num_nans>0
        % Randomly sample the model_classes labels with replacement
        sub_labels = randsample(model_classes,num_nans,true);
        classification(isnan(classification)) = sub_labels;
    end
end


end