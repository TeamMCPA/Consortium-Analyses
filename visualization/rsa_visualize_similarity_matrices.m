function rsa_visualize_similarity_matrices(model_data, test_data, model_mat, test_mat)
%% create a visualization of the (dis)similarity matrix for each subject and each session
% inputs: - model data: feature matrix for training data (condition x
%                       feature x session x subject)
%         - test_data: feature matrix for test data (condition x feature x
%                      session x subject)
%         - model_mat: similarity matrices for model data (condition x
%                      condition x session x subject)
%         - test_mat: similarity matrices for test data (condition x
%                      condition x session x subject)
% returns: a figure containing the (dis)similarity matrix for each subject and
%          each session

    % plot the training data
    plot_idx = 1;
    figure
 
    for session_idx = 1:size(model_data,3)        
        for subject_idx = 15%1:size(model_data,4)
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
    
    % plot the test data
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
