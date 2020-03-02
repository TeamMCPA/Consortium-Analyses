%% Pattern manipulation testing script
% This script builds an MCPA struct from real MCP data and then replaces
% the pattern matrix with stand-in values that are easy to visually inspect
% (e.g., integers 0-9). Then the script steps through different
% concatenation and summarization operations to allow visual confirmation
% that these functions are behaving as expected on the pattern matrix.

%% Prep workspace for the analysis
% Change to working directory where .mcp files are stored
cd('~/Documents/Projects/Consortium/Pton_adult_repairNames');

% Import the MCP file and clean up duplicate trials (a bug specific to the
% present experiment, not generally necessary)
load('SemanticNIRS_Adult_dataset1_21-Nov-2019.mcp','-mat');
MCP_data = clean_shimadzu_trials(MCP_data);

% Load the Consortium tools into the workspace
addpath(genpath('~/Documents/MATLAB/repos/Consortium-Analyses'));

%% PREVIEW OF THE DATA USING MCPA FUNCTION
% This section is only for plotting the mean HRFs and visualizing the
% overall similarity of the stimulus classes. It's not necessary to do this
% for classification. (It happens automatically in the analyses below.)

% Extract events into MCPA struct and graph the grand-average HRF over time
% for each condition. Then summarize the data by taking the time-window
% average for each condition, and plot a similarity matrix

MCPA_data = MCP_to_MCPA(...
    MCP_data,...                % MCP data struct
    [1:3],...                      % Subjects to include (blank indicates all)
    [1:5],...                      % Channels to include (blank indicates all)
    [-5,30],...                 % Analysis time window specified in seconds
    [-5,0]...                   % Baseline time window specified in seconds
    );

%% Start reassigning values in the pattern matrix

pattern_dims = size(MCPA_data.patterns);

% Time model: We'll use a simple sine wave to confirm that temporal
% relationships are preserved (no discontinuities in the wave) and so that
% the average (after running summarize) is 0+c where c is a constant
% integer that we'll induce later.
time_dim = find(strcmp('time',MCPA_data.dimensions));
time_model = sin([1:pattern_dims(time_dim)]/10)';
plot([1:size(MCPA_data.patterns,time_dim)],time_model)
title('Time Model to be used')

% Impose time model into the patterns
MCPA_data.patterns = repmat(time_model,[1,pattern_dims(2:end)]);

% Repetition: Repetition will be coded as the mean of the data in 0.1s
rep_dim = find(strcmp('repetition',MCPA_data.dimensions));
for rep = 1:pattern_dims(rep_dim)
    MCPA_data.patterns(:,:,:,rep,:,:) = ...
        (MCPA_data.patterns(:,:,:,rep,:,:) - mean(time_model)) + (rep*0.1);
end

% Session: Session will be coded as the mean of the data in 10s
sess_dim = find(strcmp('session',MCPA_data.dimensions));
for sess = 1:pattern_dims(sess_dim)
    MCPA_data.patterns(:,:,:,:,sess,:) = ...
        (MCPA_data.patterns(:,:,:,:,sess,:) - mean(time_model)) + (sess*10);
end

plot(squeeze(MCPA_data.patterns(:,1,1,1,[1:pattern_dims(sess_dim)],1)))
hold on
plot(squeeze(MCPA_data.patterns(:,1,1,[1:pattern_dims(rep_dim)],1,1)))
hold off

%% Summarization & Concatenation for most classifiers

%% Between subjects
summarize_dimensions = {'repetitionXsession', 'time'};
final_dimensions = {'conditionXrepetition+sessionXsubject', 'feature'};
summary_handle = @mean;

MCPA_summ = summarize_MCPA_Struct(summary_handle,...
    MCPA_data,...
    summarize_dimensions);
MCPA_summ.dimensions
size(MCPA_summ.patterns)

figure;
subplot(1,3,1)
plot(squeeze(MCPA_summ.patterns(1,1,:,1)),'o')
title('After summarization (1 subject)')

[group_data, group_labels, subj_data, subj_labels] = split_test_and_train(1,...
    MCPA_summ.event_types,...
    MCPA_summ.patterns,...
    MCPA_summ.event_types,...
    final_dimensions,...
    MCPA_summ.dimensions, [], []);
final_dimensions
size(subj_data)
size(group_data)

subplot(1,3,2)
plot(subj_data(:,1),'o')
title('After concatenation (test subject)')

subplot(1,3,3)
plot(group_data(:,1),'o')
title('After concatenation (training subjects)')

%% Within Subjects
summarize_dimensions = {'time'};
final_dimensions = {'conditionXrepetitionXsession', 'feature'};
summary_handle = @mean;

MCPA_summ = summarize_MCPA_Struct(summary_handle,...
    MCPA_data,...
    summarize_dimensions);
MCPA_summ.dimensions
size(MCPA_summ.patterns)

figure;
subplot(1,3,1)
plot(squeeze(MCPA_summ.patterns(1,1,:,:,1)),'o')
title('After summarization (1 subject)')

[group_data, group_labels, subj_data, subj_labels] = split_test_and_train(1,...
    MCPA_summ.event_types,...
    MCPA_summ.patterns,...
    MCPA_summ.event_types,...
    final_dimensions,...
    MCPA_summ.dimensions, [], []);
final_dimensions
size(subj_data)
size(group_data)

subplot(1,3,2)
plot(subj_data(:,1),'o')
title('After concatenation (test session)')

subplot(1,3,3)
plot(group_data(:,1),'o')
title('After concatenation (training sessions)')

%% rsa_classify Between Subjects
summarize_dimensions = {'repetition', 'time'};
final_dimensions = {'condition', 'feature', 'sessionXsubject'};
summary_handle = @mean;

MCPA_summ = summarize_MCPA_Struct(summary_handle,...
    MCPA_data,...
    summarize_dimensions);
MCPA_summ.dimensions
size(MCPA_summ.patterns)

figure;
subplot(1,3,1)
plot(squeeze(MCPA_summ.patterns(1,1,:,1)),'o')
title('After summarization (1 subject)')

[group_data, group_labels, subj_data, subj_labels] = split_test_and_train(1,...
    MCPA_summ.event_types,...
    MCPA_summ.patterns,...
    MCPA_summ.event_types,...
    final_dimensions,...
    MCPA_summ.dimensions, [], []);
final_dimensions
size(subj_data)
size(group_data)

subplot(1,3,2)
plot(squeeze(subj_data(1,1,:)),'o')
title('After concatenation (test subject)')

subplot(1,3,3)
plot(squeeze(group_data(1,1,:)),'o')
title('After concatenation (training subjects)')

%% rsa_classify Within Subjects
summarize_dimensions = {'repetition', 'time'}; 
        final_dimensions = {'condition', 'feature', 'session'};
summary_handle = @mean;

MCPA_summ = summarize_MCPA_Struct(summary_handle,...
    MCPA_data,...
    summarize_dimensions);
MCPA_summ.dimensions
size(MCPA_summ.patterns)

figure;
subplot(1,3,1)
plot(squeeze(MCPA_summ.patterns(1,1,:,1)),'o')
title('After summarization (1 subject)')

[group_data, group_labels, subj_data, subj_labels] = split_test_and_train(1,...
    MCPA_summ.event_types,...
    MCPA_summ.patterns,...
    MCPA_summ.event_types,...
    final_dimensions,...
    MCPA_summ.dimensions, [], []);
final_dimensions
size(subj_data)
size(group_data)

subplot(1,3,2)
plot(squeeze(subj_data(1,1,:)),'o')
title('After concatenation (test subject)')

subplot(1,3,3)
plot(squeeze(group_data(1,1,:)),'o')
title('After concatenation (training subjects)')

%% mcpa_classify Between subjects

% Between subjects
summarize_dimensions = {'repetitionXsession', 'repetition+session', 'time'};
final_dimensions = {'conditionXsubject', 'feature'};
summary_handle = @mean;

MCPA_summ = summarize_MCPA_Struct(summary_handle,...
    MCPA_data,...
    summarize_dimensions);
MCPA_summ.dimensions
size(MCPA_summ.patterns)

figure;
subplot(1,3,1)
plot(squeeze(MCPA_summ.patterns(:,1,1)),'o')
title('After summarization (1 subject)')

[group_data, group_labels, subj_data, subj_labels] = split_test_and_train(1,...
    MCPA_summ.event_types,...
    MCPA_summ.patterns,...
    MCPA_summ.event_types,...
    final_dimensions,...
    MCPA_summ.dimensions, [], []);
final_dimensions
size(subj_data)
size(group_data)

subplot(1,3,2)
plot(squeeze(subj_data(:,1)),'o')
title('After concatenation (test subject)')

subplot(1,3,3)
plot(squeeze(group_data(:,1)),'o')
title('After concatenation (training subjects)')

%% mcpa_classify Within subjects

% Between subjects
summarize_dimensions = {'repetition', 'time'}; 
        final_dimensions = {'conditionXsessionXsubject', 'feature'};
summary_handle = @mean;

MCPA_summ = summarize_MCPA_Struct(summary_handle,...
    MCPA_data,...
    summarize_dimensions);
MCPA_summ.dimensions
size(MCPA_summ.patterns)

figure;
subplot(1,3,1)
plot(squeeze(MCPA_summ.patterns(1,1,:,1)),'o')
title('After summarization (1 subject)')

[group_data, group_labels, subj_data, subj_labels] = split_test_and_train(1,...
    MCPA_summ.event_types,...
    MCPA_summ.patterns,...
    MCPA_summ.event_types,...
    final_dimensions,...
    MCPA_summ.dimensions, [], []);
final_dimensions
size(subj_data)
size(group_data)

subplot(1,3,2)
plot(squeeze(subj_data(:,1)),'o')
title('After concatenation (test subject)')

subplot(1,3,3)
plot(squeeze(group_data(:,1)),'o')
title('After concatenation (training subjects)')


