function HRF_builder( subj_list )

%% Get HRFs for subject
% first extract_vecs_1subj - get windowed data for each subject
    % this includes another function: scrub_artifacts

for i = 1:length(subj_list)
    try
    seek_files = dir([subj_list{i} '/' subj_list{i} '*.nirs']);
    subjs(i).Condition = extract_vecs_1subj(...
        [seek_files(1).folder filesep seek_files(1).name],...
        {'ball','cookie','cup','shoe'},[0,165],'none');
    catch
        disp(['CANNOT FIND SUBJECT ' subj_list{i}]);
    end
end