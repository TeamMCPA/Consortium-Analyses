function classification = svm_classify(train_data, train_labels, test_data, test_labels, opts)
%% svm_classify implements support vector machine 
% input: model_data - data from the subjects in the group model
%        test_data -  labels for the subjects in the group model
%        test_data - subject data we want to predict for
%        test_labels - subject labels, unused
%        opts (optional) - an options struct

%% create SVM model based on opts struct

% check to see if we need to call fitcecoc or fitcsvm
if isfield(opts, 'Multiclass')
    multiclass = opts.Multiclass;
    opts = rmfield(opts,'Multiclass');
else
    multiclass = false;
end

% parse training parameters
input = parse_opts(opts);

% build model
if multiclass
    t = templateSVM(input{:});
    svm_model = fitcecoc(train_data, train_labels,'Learners', t);
else
    svm_model = fitcsvm(train_data, train_labels, input{:});
end



% classify the participant
classification = predict(svm_model, test_data);

end