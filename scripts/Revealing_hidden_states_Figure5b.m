clear all;
close all;
clc;
%% Script for figure 5b "Revealing hidden states in visual working memory using electroencephalography"

root_dir =  'your path'; % path to data and functions 
addpath(genpath(root_dir));

angle_bins= [1 2 3 4];
cond_pairs = [1 3;2 4];
test_chans = { 'P7','P5','P3','P1','Pz','P2','P4','P6','P8','PO7','PO3','POz','PO4',...
    'PO8','O2','O1','Oz'}'; % channels to be used in analysis
%% Mahalanobis Distance of the
for sub=1:24
    fprintf(['Doing ' num2str(sub) '\n'])
    if sub<10
        load(fullfile(root_dir,['Data_WM_hidden_states_s0' num2str(sub) '.mat']));
    else
        load(fullfile(root_dir,['Data_WM_hidden_states_s' num2str(sub) '.mat']));
    end
    % exclude bad trials
    incl_imp=setdiff(1:length(Long_trials.EEG_impulse_only_tm_impulse.trial),Long_trials.bad_trials);
    
    % extract only relevant time-points, channels, smooth data, and put in
    % the right format
    t=1; c=1;
    Impulse=zeros(length(incl_imp),length(test_chans),length(find(Long_trials.EEG_impulse_only_tm_impulse.time{1}<=0.8)));
    for trl=incl_imp
        for ch=find(ismember(Long_trials.EEG_impulse_only_tm_impulse.label,test_chans))
            Impulse(t,c,:)=gsmooth(Long_trials.EEG_impulse_only_tm_impulse.trial{trl}(ch,Long_trials.EEG_impulse_only_tm_impulse.time{1}<=0.8),2);
            c=c+1;
        end
        t=t+1;c=1;
    end
    
    % compute trial wise mahalanobis distance difference between orthogonal angle bins
    for pair=1:size(cond_pairs,1)
        dat1=Impulse(ismember(Long_trials.conditions.item_angle_bins1(incl_imp,1),angle_bins(cond_pairs(pair,1)))==1,:,:);
        dat2=Impulse(ismember(Long_trials.conditions.item_angle_bins1(incl_imp,1),angle_bins(cond_pairs(pair,2)))==1,:,:);
        d = mahalFuncTrialWise([dat1;dat2],[ones(size(dat1,1),1);ones(size(dat2,1),1)*2]);
        d1_p_imp(pair,:)=mean(d,1);
    end
    for pair=1:size(cond_pairs,1)
        dat1=Impulse(ismember(Long_trials.conditions.item_angle_bins2(incl_imp,1),angle_bins(cond_pairs(pair,1)))==1,:,:);
        dat2=Impulse(ismember(Long_trials.conditions.item_angle_bins2(incl_imp,1),angle_bins(cond_pairs(pair,2)))==1,:,:);
        d = mahalFuncTrialWise([dat1;dat2],[ones(size(dat1,1),1);ones(size(dat2,1),1)*2]);
        d2_p_imp(pair,:)=mean(d,1);
    end
    Impulse_mahal(:,sub)=mean([d1_p_imp;d2_p_imp],1)'; %average across all trials
end
time=Long_trials.EEG_impulse_only_tm_impulse.time{1}(Long_trials.EEG_impulse_only_tm_impulse.time{1}<=0.8);

% test for significant decoding after impulse presentation
[datobs_impulse, datrnd_impulse] = cluster_test_helper(Impulse_mahal, 10000,'t');
[h_impulse, p_impulse, ~] = cluster_test(datobs_impulse(time>=0,1), datrnd_impulse(time>=0,:), 0, 0.05,0.01);

% test for significant decoding increase after impulse presentation
Impulse_mahal_increase=Impulse_mahal-repmat(mean(Impulse_mahal(time>=-0.1 & time<=0,:),1),[size(Impulse_mahal,1),1]);
[datobs_impulse_incr, datrnd_impulse_incr] = cluster_test_helper(Impulse_mahal_increase, 10000,'t');
[h_impulse_incr, p_impulse_incr, ~] = cluster_test(datobs_impulse_incr(time>=0,1), datrnd_impulse_incr(time>=0,:), 0, 0.05,0.01);

figure
line('XData', [min(time) max(time)], 'YData', [0 0], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color','k');
hold all
plot(time,mean(Impulse_mahal,2))
scatter(time([false(1,length(time(time<0))) h_impulse']),ones(1,length(find(h_impulse==1))).*0.0005,'b','filled')
scatter(time([false(1,length(time(time<0))) h_impulse_incr']),ones(1,length(find(h_impulse_incr==1))).*0.012,'MarkerEdgeColor',[0 0 0.5],'MarkerFaceColor',[0 0 0.5])
ylim([-0.002 0.04]);
pbaspect([1,1.25,1])
xlabel('Time relative to impulse')
ylabel('Multivaraite discrimination')

