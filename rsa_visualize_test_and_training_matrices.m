function rsa_visualize_test_and_training_matrices(training_matrix, model_labels, test_matrix, test_labels)
%% visualize the final test and training models
% inputs: - training_matrix: (dis)similarity matrix for training model
%         - model_labels: training model labels
%         - test_matrix: (dis)similarity matrix for test model
%         - test_labels: test model labels
% returns: visualization of the test and train models


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
