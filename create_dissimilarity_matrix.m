function dissimilarity_matrix = create_dissimilarity_matrix(feature_matrix, opts)
    dissimilarity_matrix = nan(size(feature_matrix,1),size(feature_matrix,1),size(feature_matrix,3),size(feature_matrix,4));
    for i = 1: (size(feature_matrix,3)*size(feature_matrix,4))
        dissimilarity_matrix(:,:,i) = squareform(pdist(feature_matrix(:,:,i), opts.metric));
    end
end