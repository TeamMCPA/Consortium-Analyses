function allsubj_results = leave_one_Ss_out_classifyIndividualEvents2(MCPA_pattern,cond1,cond2,setsize)

%% leave_one_Ss_out_classifyIndividualEvents2 takes a struct output from the 
% estimate_WindowAverages function and performs the leave-one-subject-out
% test for each subject to classify every event independently for the
% left-out (test) subject. Multichannel patterns for each event are
% compared to the group model and classified accordingly.
%
% The function is called with the following parameters:
% leave_one_Ss_out_classifyIndividualEvents( MCPA_pattern , cond1 , cond2 [, setsize])
%
% MCPA_pattern: the output of the summarize_MCPA_Struct function. This
% structure should contain the metadata about how multichannel patterns
% were estimated (which subjects and channels were included) and the
% multichannel patterns themselves. 
%
% cond1 and cond2: the event types to be compared (corresponding to the
% rows in MCPA_pattern.pattern, not necessarily the original mark number
% from the HomER file). If there are only 2 event types, cond1 and cond2
% may simply be 1 and 2, respectively. If there are 4 event types (e.g.,
% cat, dog, shoe, boot), then single conditions maybe be selected as a 
% subset (e.g., cond1=3 and cond2=4 for shoe vs. boot) or as vectors (e.g.,
% cond1=[1,2] and cond2=[3,4] for animals vs. footwear).
%
% setsize (optional): the number of channels that will be included in the
% analysis. If setsize is less than the total number of channels in
% MCPA_pattern struct (MCPA_pattern.incl_channels), then all combinations
% of the channels will be analized in a n-choose-k fashion. In this case,
% average decoding accuracy across combinations will be computed for each
% channel. Default behavior is to include all channels.

if ndims(MCPA_pattern.patterns)==4
    handle = input('MCPA pattern has not been summarized. Enter function name to apply for summarization: ','s');
    MCPA_pattern = summarize_MCPA_Struct(handle,MCPA_pattern);
end
subj_mat = MCPA_pattern.patterns;

%% Check the setsize, if it is specified at all.
if nargin < 4
    % default for setsize is all channels
    % notice the dimension of subj_mat is (event types x channels x subjects)
    % setsize = number of channels
    setsize = size(subj_mat,2);
else
    if setsize > size(subj_mat,2)
        warning('Set size is larger than the total number of channels. Ignoring.');
        setsize = size(subj_mat,2);
    end
    fprintf('Running analysis based on a set size of %d channels...\n', setsize);
end


%% Classifying individual events
fprintf('\nLabel by Event (1x2 similarity matrix for conditions [');
fprintf(' %d', cond1);
fprintf(' ] vs. [');
fprintf(' %d', cond2);
fprintf(' ])\n');


%% Load in the specified data file as "nirs_data" and then rename the data struct as "nirs_data"
fprintf('Loading data from the data file...');

if length(MCPA_pattern.data_file) == 1  % The case that input is an MCP file.
    [~, ~, ext] = fileparts(MCPA_pattern.data_file{:});
    if strcmp(ext,'.mcp')
        mcp_file_content = load(data_file,'-mat');
        varname = fieldnames(mcp_file_content);
        nirs_data = eval(['mcp_file_content.' varname{1}]);
        clear('mcp_file_content')
    else
        error('data_file is not MCP. Will not be able to load trial-level data');
    end
    
elseif length(MCPA_pattern.data_file) > 1 % The case that input need to be converted to MCP
    nirs_data = build_MCP(MCPA_pattern.data_file,cellstr(num2str(MCPA_pattern.incl_subjects'))',repmat({'probe'},1,length(MCPA_pattern.incl_subjects)),'s');
    
else 
    error('There is something wrong with the file input...');
end

%% Setting up the combinations of channel subsets
% Create all possible subsets. If setsize is equal to the total number of
% channels, there will only be one 'subset' which is the full channel
% array. If setsize is less than the total number of channels, there will
% be n-choose-k subsets to analyze.
sets = nchoosek(1:size(subj_mat,2),setsize);
num_subjs = size(subj_mat,3);

%% Set up the results structure which includes a copy of MCPA_pattern

allsubj_results = [];
% save some basic info about in the allsubj_results
allsubj_results.MCPA_pattern = MCPA_pattern;
allsubj_results.data_file = MCPA_pattern.data_file;
allsubj_results.created = datestr(now);
allsubj_results.test_type = 'Leave one subject out, Classify individual events';
allsubj_results.setsize = setsize;
allsubj_results.incl_subjects = MCPA_pattern.incl_subjects;
allsubj_results.incl_channels = MCPA_pattern.incl_channels;
allsubj_results.cond1 = cond1;
allsubj_results.cond2 = cond2;
allsubj_results.subsets = MCPA_pattern.incl_channels(sets);

allsubj_results.accuracy.cond1.subjXchan = nan(num_subjs,size(subj_mat,2));
allsubj_results.accuracy.cond2.subjXchan = nan(num_subjs,size(subj_mat,2));
allsubj_results.accuracy.cond1.subsetXsubj = nan(num_subjs,size(sets,1))';
allsubj_results.accuracy.cond2.subsetXsubj = nan(num_subjs,size(sets,1))';

%% Iterate through all the subjects and attempt classification

if setsize==size(subj_mat,2)
    fprintf('  Subject\tCond1\tCond2\n');
else
    fprintf('  Subject-wise results suppressed for subset analysis.\n');
end

for s_idx = 1:num_subjs
    
    % Group holds the _indexes_ of the included subjects, not their actual
    % subject numbers. Drop the index for the subject who is going to be
    % tested in this iteration of the loop (i is his/her index)
    group = [1:num_subjs];
    group(s_idx) = [];
    
    % Build the averaged data for the group (non-test) subjects so that the
    % first dimension is cond1 and cond2 accordingly
    group_vecs = mean(subj_mat(:,:,group),3)';
    group_vecs = cat(2,mean(group_vecs(:,cond1),2),mean(group_vecs(:,cond2),2));
    
    %% Extract oxy data from the MCP struct for test subject
    events = MCP_get_subject_events(nirs_data(s_idx),1:size(subj_mat,2),MCPA_pattern.time_window,{['m' num2str(cond1)],['m' num2str(cond2)]});

    for currentmark = 1:2 % Only two conditions allowed at the moment
        
        % to store the results for the number of channels, by the number of
        % channels, by whether or not the channel was accurate for that search
        % light as well as the accurcies for the matching vs. the mismatching
        % conditions
        results(currentmark).acc = nan(size(sets,1),size(subj_mat,2),2);
        results(currentmark).markNo = currentmark;
        
        for search = 1:size(sets,1)
            
            for k = 1:size(events,3)
                % Only analyze if there are data for this repeition k of
                % the currentmark
                if sum(sum(~isnan(events(:,:,k,currentmark))))>0                 
                    if setsize > 1
                        % group_vecs(:,1) is always for cond1, group_vecs(:,2) is
                        % always for cond2
                        gtype1 = corr(group_vecs(sets(search,:),1),nanmean(events(:,sets(search,:),k,currentmark),1)');
                        gtype2 = corr(group_vecs(sets(search,:),2),nanmean(events(:,sets(search,:),k,currentmark),1)');
                        
                        % classify the accuracy based on which condition this mark
                        % is part of
                        if currentmark==1
                            results(currentmark).acc(search,sets(search,:),k) = gtype1>gtype2;
                        elseif currentmark==2
                            results(currentmark).acc(search,sets(search,:),k) = gtype1<gtype2;
                        end
                    elseif setsize == 1
                        gtype1 = abs(group_vecs(sets(search,:),1)-mean(events(:,sets(search,:),k,currentmark),1)');
                        gtype2 = abs(group_vecs(sets(search,:),2)-mean(events(:,sets(search,:),k,currentmark),1)');
                        
                        % classify the accuracy based on which condition this mark
                        % is part of, in the case where the setsize is 1 we
                        % are looking for the closest value not the largest
                        % correlation
                        if sum(totalmarks(currentmark)==cond1marks)>0
                            results(currentmark).acc(search,sets(search,:),k) = gtype1<gtype2;
                        elseif sum(totalmarks(currentmark)==cond2marks)>0
                            results(currentmark).acc(search,sets(search,:),k) = gtype1>gtype2;
                        end
                    end
                end
            end
        end
    end
    
    % determine where to put this accuracy data based on the
    % conditions, i.e. matching the marks to the conditions
    
    acc_cond1 = [];
    acc_cond2 = [];
    
    for mark = 1:length(totalmarks)
        if sum(totalmarks(mark)==cond1marks)>0
            acc_cond1 = cat(1,acc_cond1,events(mark).acc);
        elseif sum(totalmarks(mark)==cond2marks)>0
            acc_cond2 = cat(1,acc_cond2,events(mark).acc);
        end
    end
    
    % calculate the mean for all sets for this subject before the data is
    % lost forever! (*quieter now, like a fading echo* forever, forever, forever, forever...)
    setresults = nan(size(sets,1),1);
    foo = vertcat(events(1).acc);
    for ii = 1:size(sets,1)
        setresults(ii) = nanmean(foo(:,sets(ii,1),ii));
    end
    allsubj_results.accuracy.cond1.subsetXsubj(:,s_idx)=setresults';
    
    setresults = nan(size(sets,1),1);
    foo = vertcat(events(2).acc);
    for ii = 1:size(sets,1)
        setresults(ii) = nanmean(foo(:,sets(ii,1),ii));
    end
    allsubj_results.accuracy.cond2.subsetXsubj(:,s_idx)=setresults';

    % if we're just running this thing once.. just spit out the answer
    if setsize == size(subj_mat,2)
        fprintf('  %d\t',incl_subjects(s_idx));
        fprintf('\t%2.0f%%\t%2.0f%%\n',nanmean(acc_cond1(:,1))*100,nanmean(acc_cond2(:,1))*100);
        allsubj_results.accuracy.cond1.subjXchan(s_idx,:) = nanmean(acc_cond1,1);
        allsubj_results.accuracy.cond2.subjXchan(s_idx,:) = nanmean(acc_cond2,1);
    else
        allsubj_results.accuracy.cond1.subjXchan(s_idx,:) = nanmean(nanmean(acc_cond1,3),1);
        allsubj_results.accuracy.cond2.subjXchan(s_idx,:) = nanmean(nanmean(acc_cond2,3),1);
    end
    
    
    
end

if setsize ~= size(subj_mat,2)

    figure
    errorbar(1:size(allsubj_results.accuracy.cond1.subjXchan,2),mean(allsubj_results.accuracy.cond1.subjXchan),std(allsubj_results.accuracy.cond1.subjXchan)/sqrt(size(allsubj_results.accuracy.cond1.subjXchan,1)),'r')
    hold;
    errorbar(1:size(allsubj_results.accuracy.cond2.subjXchan,2),mean(allsubj_results.accuracy.cond2.subjXchan),std(allsubj_results.accuracy.cond2.subjXchan)/sqrt(size(allsubj_results.accuracy.cond2.subjXchan,1)),'k')
    title('Decoding Accuracy across all channels: Red = Cond1, Black = Cond2')
    set(gca,'XTick',[1:length(incl_channels)])
    set(gca,'XTickLabel',incl_channels)
    hold off;
    
    figure
    errorbar(1:size(allsubj_results.accuracy.cond1.subjXchan,1),mean(allsubj_results.accuracy.cond1.subjXchan'),repmat(std(mean(allsubj_results.accuracy.cond1.subjXchan'))/sqrt(size(allsubj_results.accuracy.cond1.subjXchan,2)),1,size(allsubj_results.accuracy.cond1.subjXchan,1)),'r')
    hold;
    errorbar(1:size(allsubj_results.accuracy.cond2.subjXchan,1),mean(allsubj_results.accuracy.cond2.subjXchan'),repmat(std(mean(allsubj_results.accuracy.cond2.subjXchan'))/sqrt(size(allsubj_results.accuracy.cond2.subjXchan,2)),1,size(allsubj_results.accuracy.cond2.subjXchan,1)),'k')
    title('Decoding Accuracy across all subjects: Red = Cond1, Black = Cond2')
    set(gca,'XTick',[1:length(incl_subjects)])
    set(gca,'XTickLabel',incl_subjects)
    hold off;
end
