function classification = logit_classify(train_data, train_labels, test_data, test_label, opts)
%% logistic regression wrapper
% takes in training data, training labels, testing data, and testing labels
% (unused) as well as opts struct with classification parameters

% parse parameters
input = parse_opts(opts);

% create logistic regressiontemplate
template = templateLinear('Learner', 'logistic', input{:});
% build the model
model = fitcecoc(train_data, train_labels);
% classify
classification = predict(model, test_data);

end