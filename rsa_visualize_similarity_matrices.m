function rsa_visualize_similarity_matrices(model_data, test_data, model_mat, test_mat)
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