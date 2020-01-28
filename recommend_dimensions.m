function [summarize_dimensions, final_dimensions] = recommend_dimensions(results_struct, isWithinSubjects)
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
        final_dimensions = {'conditionXsessionXsubject', 'channel'};
    else
        summarize_dimensions = {'repetitionXsession', 'repetition+session', 'time'};
        final_dimensions = {'conditionXsubject', 'channel'};
    end    
elseif strcmp(func2str(results_struct.test_handle),'rsa_classify')
    if isWithinSubjects
        summarize_dimensions = {'repetition', 'time'}; 
        final_dimensions = {'condition', 'channel', 'session'};
    else
        summarize_dimensions = {'repetition', 'time'}; 
        final_dimensions = {'condition', 'channel', 'sessionXsubject'};
    end
else
    if isWithinSubjects
        summarize_dimensions = {'time'};
        final_dimensions = {'conditionXrepetitionXsession', 'channel'};
    else
        summarize_dimensions = {'repetitionXsession', 'time'};
        final_dimensions = {'conditionXrepetition+sessionXsubject', 'channel'};
    end     
end
end