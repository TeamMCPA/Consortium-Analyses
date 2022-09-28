function draw_mcpa_output( between_subj_level )
%% draw_mcpa_output visualizes output from classification scripts
% The only input to the draw_mcpa_output function is the struct resulting
% from a classification, e.g,. output from nfold_classify_ParticipantLevel.
% The function requires this struct to have fields for event_types, 
% conditions, patterns (or summed_mcpa_patterns), hemoglobin, and
% accuracy_matrix. It will output a figure illustrating the average feature
% values for each stimulus class and another figure with the pairwise
% classification accuracies for all stimulus classes.
% bdz - summer 2021

%% Plot responses patterns
% Copy the response patterns from the MCPA output
[conditions, cond_idx, event_idx] = intersect(between_subj_level.conditions, between_subj_level.event_types, 'stable');
num_cond = length(between_subj_level.conditions);

% Get information about the multivariate pattern data used in the classifier
try
    % Updated field names August 2021
    num_subj = size(between_subj_level.patterns,find(strcmp(between_subj_level.dimensions,'subject')));
    num_feat = size(between_subj_level.patterns,find(strcmp(between_subj_level.dimensions,'feature'))); 
catch
    % Old version of toolbox had different field names
    num_subj = size(between_subj_level.summed_mcpa_patterns,find(strcmp(between_subj_level.dimensions,'subject')));
    num_feat = size(between_subj_level.summed_mcpa_patterns,find(strcmp(between_subj_level.dimensions,'feature')));
end
 
% Pare down to only oxygenated hemoglobin
try
    % Updated field name 'pattern' August 2021
    hemo_patterns = between_subj_level.patterns(...
        event_idx(cond_idx),...
        1 : end,...
        :);
catch
    % Old version of toolbox had different field name 'summed_mcpa_pattern'
    hemo_patterns = between_subj_level.summed_mcpa_patterns(...
        event_idx(cond_idx),...
        1 : end,...
        :);
end

% Drop parcels where fewer than 50% of participants have data
% (for visualization purposes only)
hemo_patterns_chans = find(mean(isnan(hemo_patterns(1,:,:)),3)<.5);
oxy_patterns_filt = hemo_patterns(:,hemo_patterns_chans,:);
% Transform into a 2d matrix by condition and then by subject
oxy_patterns_sorted = reshape(permute(oxy_patterns_filt,[1 3 2]), num_cond*num_subj, size(oxy_patterns_filt,2));
oxy_patterns_labels = repmat(between_subj_level.conditions',num_subj,1);

% Average response patterns
figure;
imagesc(nanmean(oxy_patterns_filt,3))
if size(oxy_patterns_filt,2)<50, axis image; end;
xlabel('Feature Number');
% xticks([1:length(hemo_patterns_chans)]);
% xticklabels(hemo_patterns_chans);
yticks([1:num_cond]);
yticklabels(between_subj_level.conditions);
title(['Average Response Patterns: ' between_subj_level.hemoglobin 'genated Hemoglobin'])
if strcmpi(between_subj_level.hemoglobin,'Deoxy')
    colormap('cool');
elseif strcmpi(between_subj_level.hemoglobin,'Oxy')
    colormap('hot');
else
    colormap('parula');
end
colorbar('eastoutside');
fig = gca;
fig.FontSize = 16;

% Average Classification accuracy
if length(between_subj_level.conditions)>2
    figure;
    acc = squeeze(nanmean(between_subj_level.accuracy_matrix,4))';
    imagesc( acc )
    xticks(1:num_cond);
    xticklabels(between_subj_level.conditions);
    yticks(1:num_cond);
    yticklabels(between_subj_level.conditions);
    caxis([0,1]);
    colorbar('eastoutside');
    colormap('parula');
    el_id = 1;
    for x=1:num_cond
        for y=1:num_cond
            text(x,y,num2str(round(acc(el_id),2)),'FontSize',16);
            el_id = el_id+1;
        end
    end   
    fig = gca;
    fig.FontSize = 16;
end

% % MDS Plot of subject-level patterns
% figure;
% oxy_patterns_sorted(isnan(oxy_patterns_sorted)) = 0;
% D = pdist(oxy_patterns_sorted,'spearman');
% [Y e] = cmdscale(D);
% markcol = {'r','b','g','k'};
% for i = 1:num_cond
%     marker = [markcol{i} 'o'];
%     plot(Y(i:num_cond:end,1),Y(i:num_cond:end,2),marker,'MarkerSize',10,'LineWidth',3)
%     hold on;
% end
% legend(between_subj_level.conditions)
% title('2D Subject-Level Peekaboo Response Patterns','FontSize',14)
% fig = gca;
% fig.FontSize = 16;
% hold off;
