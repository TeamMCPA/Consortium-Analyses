function similarity_matrix = create_similarity_matrix(feature_matrix, opts)
%% abstract feature vectors into similarity space
% inputs: -feature matrix (feature representations with dimensions:
%         condition x feature x session x subject
%         -opts (your opts struct)
% returns: matrix containing the similarity structure for each session
%          and each subject
similarity_matrix = nan(size(feature_matrix,1),size(feature_matrix,1),size(feature_matrix,3),size(feature_matrix,4));
for i = 1: (size(feature_matrix,3)*size(feature_matrix,4))
    % want to remove completely empty columns, but not empty rows (to
    % preserve condition order)
    keep_inds = sum(ismissing(feature_matrix(:,:,i))) ~= size(feature_matrix(:,:,i),1);
    if sum(keep_inds) > 0
        feature_mat = feature_matrix(:,keep_inds,i);
    else
        feature_mat = feature_matrix(:,:,i);
    end

    similarity_matrix(:,:,i) = corr(feature_mat', 'type', opts.metric);

end
end