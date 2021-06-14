function dissimilarity_matrix = create_dissimilarity_matrix(feature_matrix, opts)
%% abstract feature vectors into dissimilarity space
% inputs: -feature matrix (feature representations with dimensions:
%         condition x feature x session x subject
%         -opts (your opts struct)
% returns: matrix containing the dissimilarity structure for each session
%          and each subject

    dissimilarity_matrix = nan(size(feature_matrix,1),size(feature_matrix,1),size(feature_matrix,3),size(feature_matrix,4));
    for i = 1: (size(feature_matrix,3)*size(feature_matrix,4))
        dissimilarity_matrix(:,:,i) = squareform(pdist(feature_matrix(:,:,i), opts.metric));
    end
end
