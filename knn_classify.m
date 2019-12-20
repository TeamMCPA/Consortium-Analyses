function classification = knn_classify(train_data, train_labels, test_data, test_labels, opts)
%% wrapper for KNN classification
% takes in training data, training labels, testing data, testing labels and
% an opts struct with parameters for classification
% outputs the classification results

%% parse out the classification parameters
if isstruct(opts) 
    input = parse_opts(opts);
    knn_model = fitcknn(train_data, train_labels, input{:});
else
    knn_model = fitcknn(train_data, train_labels);

end

%% predict
classification = predict(knn_model, test_data);

end