function summarize_dimensions = recommend_summarize_dimensions(results_struct, isWithinSubjects)
%% this function recommends dimensions for different methods of classification
% input:
% results_struct - input into our cross validation function. We use this to
% determine what kind of classification the user wants
% isWithinSubjects - marks whether the user is doing within subjects cross
% validation. This would call for different recommended dimension than
% between subjects cross validation

% output: a list of the recommended dimensions to summarize over

if strcmp(func2str(results_struct.test_handle),'mcpa_classify')
    if isWithinSubjects
        summarize_dimensions = {'repetition', 'time'}; 
    else
        summarize_dimensions = {'repetitionXsession', 'repetition+session', 'time'};
    end    
elseif strcmp(func2str(results_struct.test_handle),'rsa_classify')
    summarize_dimensions = {'repetition', 'time'};      
else
    if isWithinSubjects
        summarize_dimensions = {'time'};
    else
        summarize_dimensions = {'repetitionXsession', 'time'};
    end     
end
end