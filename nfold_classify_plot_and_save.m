function nfold_classify_plot_and_save(cond_names,nfold_category_results)

%% function plot and save the results from nfold_classify functions

% Arguments: 
% nfold_category_results: n_fold_category_results structure from  
%       n_fold_classify_ParticipantLevel.m or  n_fold_classify_WithinSubjects.m
% cond_names: List of the names in trigger # order (1-8)
%       e.g., cond_names = {'baby','book','bottle','cat','dog','hand','shoe','spoon'};


%% create a results folder  

if ~exist('results', 'dir')
   mkdir('results');
end

%% save the results 

filename = sprintf('results/%s_results.mat',nfold_category_results.test_type);
save(filename,'nfold_category_results');

%% partiicpant level plot

if strcmp(nfold_category_results.test_type, 'nfold_classify_ParticipantLevel')
    fig1 = figure('Position', get(0, 'Screensize'));
    mean_acc = squeeze(nanmean(nfold_category_results.accuracy_matrix,4));
    imagesc(mean_acc')
    title('participant-level accuracy plot');
    xticklabels(cond_names)
    yticklabels(cond_names)
    caxis([0,1])
    colormap('parula')
    colorbar('SouthOutside')
    [i, j, ~] = find(~isnan(mean_acc));
    text(i-.3,j,num2str(round(mean_acc(~isnan(mean_acc)),2)));
    nanmean(mean_acc(:))
    
    % save the plots 
    saveas(fig1,'results/participant_level','png');
end 

%% within participant level plots 

if strcmp(nfold_category_results.test_type, 'nfold_classify_WithinSubjects')

    % Plot the mean accuracy for each pairwise comparison for each participant
    fig2 = figure('Position', get(0, 'Screensize'));
    mean_acc = squeeze(nanmean(nfold_category_results.accuracy_matrix,4));
    for i = 1:16
        subplot(4,4,i);
        imagesc(mean_acc(:,:,i)')
        title(sprintf('subject #%d',i));
        xticklabels([])
        yticklabels([])
        caxis([0,1])
        colormap('parula')
        colorbar('SouthOutside')
    end 

    % Plot the mean accuracy for each pairwise comparison for all participants
    fig3 = figure('Position', get(0, 'Screensize'));
    mean_acc = nanmean(mean_acc,3); % average across all 16 participants
    imagesc(mean_acc')
    title('averaged accuracy across all subjects');
    xticklabels(cond_names)
    yticklabels(cond_names)
    caxis([0,1])
    colormap('parula')
    colorbar('SouthOutside')
    [i, j, ~] = find(~isnan(mean_acc));
    text(i-.3,j,num2str(round(mean_acc(~isnan(mean_acc)),2)));
    nfold_category_results.mean_accuracy = nanmean(mean_acc(:))

    % save the plots 
    saveas(fig2,'results/within_participant_plot_all_subs','png');
    saveas(fig3,'results/within_participant_plot_each_sub','png');
end 


end

