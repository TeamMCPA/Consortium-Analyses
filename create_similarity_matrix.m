function similarity_matrix = create_similarity_matrix(feature_matrix, opts)
%% abstract feature vectors into similarity space
% inputs: -feature matrix (feature representations with dimensions:
%         condition x feature x session x subject
%         -opts (your opts struct)
% returns: matrix containing the similarity structure for each session
%          and each subject
    
    similarity_matrix = nan(size(feature_matrix,1),size(feature_matrix,1),size(feature_matrix,3),size(feature_matrix,4));
    for i = 1: (size(feature_matrix,3)*size(feature_matrix,4))
        similarity_matrix(:,:,i) = corr(feature_matrix(:,:,i)','rows', 'pairwise', 'type', opts.metric);
    end
end
