clear all;
close all;
clc;

root_dir =  'your path'; % path to data and functions 
addpath(genpath(root_dir));

angle_bins= [1 2 3 4]; % the four angle bins
cond_pairs = [1 3;2 4]; %angle bin comparisons
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
    
    % extract only relevant time-points, channels, smooth data, and put in the right format
    incl_long=setdiff(1:length(Long_trials.EEG_mem_item_only.trial),Long_trials.bad_trials);
    t=1; c=1;
    M_long=zeros(length(incl_long),length(test_chans),length(Long_trials.EEG_mem_item_only.time{1}));
    for trl=incl_long
        for ch=find(ismember(Long_trials.EEG_mem_item_only.label,test_chans))
            M_long(t,c,:)=gsmooth(Long_trials.EEG_mem_item_only.trial{trl}(ch,:),2);
            c=c+1;
        end
        t=t+1;c=1;
    end
    incl_short=setdiff(1:length(Short_trials.EEG_mem_item.trial),Short_trials.bad_trials);
    t=1; c=1;
    M_short=zeros(length(incl_short),length(test_chans),length(Short_trials.EEG_mem_item.time{1}));
    for trl=incl_short
        for ch=find(ismember(Short_trials.EEG_mem_item.label,test_chans))
            M_short(t,c,:)=gsmooth(Short_trials.EEG_mem_item.trial{trl}(ch,:),2);
            c=c+1;
        end
        t=t+1;c=1;
    end
    
    % compute trial wise mahalanobis distance difference between orthogonal angle bins
    angle_bins1_long=Long_trials.conditions.item_angle_bins1(incl_long);
    angle_bins1_short=Short_trials.conditions.item_angle_bins1(incl_short);
    angle_bins2_long=Long_trials.conditions.item_angle_bins2(incl_long);
    angle_bins2_short=Short_trials.conditions.item_angle_bins2(incl_short);
    for pair=1:size(cond_pairs,1)
        dat1_long=M_long(ismember(angle_bins1_long,angle_bins(cond_pairs(pair,1)))==1,:,:);
        dat2_long=M_long(ismember(angle_bins1_long,angle_bins(cond_pairs(pair,2)))==1,:,:);
        d = mahalFuncTrialWise([dat1_long;dat2_long],[ones(size(dat1_long,1),1);ones(size(dat2_long,1),1)*2]);
        d1_p_long(pair,:)=mean(d,1);
        
        dat1_short=M_short(ismember(angle_bins1_short,angle_bins(cond_pairs(pair,1)))==1,:,:);
        dat2_short=M_short(ismember(angle_bins1_short,angle_bins(cond_pairs(pair,2)))==1,:,:);
        d = mahalFuncTrialWise([dat1_long;dat1_short;dat2_long;dat2_short],...
            [ones(size([dat1_long;dat1_short],1),1);ones(size([dat2_long;dat2_short],1),1)*2]);
        d1_p_all(pair,:)=mean(d,1);
        
        dat1_long=M_long(ismember(angle_bins2_long,angle_bins(cond_pairs(pair,1)))==1,:,:);
        dat2_long=M_long(ismember(angle_bins2_long,angle_bins(cond_pairs(pair,2)))==1,:,:);
        d = mahalFuncTrialWise([dat1_long;dat2_long],[ones(size(dat1_long,1),1);ones(size(dat2_long,1),1)*2]);
        d2_p_long(pair,:)=mean(d,1);
        
        dat1_short=M_short(ismember(angle_bins2_short,angle_bins(cond_pairs(pair,1)))==1,:,:);
        dat2_short=M_short(ismember(angle_bins2_short,angle_bins(cond_pairs(pair,2)))==1,:,:);
        d = mahalFuncTrialWise([dat1_long;dat1_short;dat2_long;dat2_short],...
            [ones(size([dat1_long;dat1_short],1),1);ones(size([dat2_long;dat2_short],1),1)*2]);
        d2_p_all(pair,:)=mean(d,1);
    end
    M_long_mahal(:,sub)=mean([d1_p_long;d2_p_long],1)';
    M_all_mahal(:,sub)=mean([d1_p_all;d2_p_all],1)';
end
time=Short_trials.EEG_mem_item.time{1};

% test for significant decoding in long trials only
[datobs_long, datrnd_long] = cluster_test_helper(M_long_mahal, 10000,'t');
[h_long, p_long, ~] = cluster_test(datobs_long(time>=0,1),datrnd_long(time>=0,:),0,0.05,0.01);

% test test for significant decoding in all trials
[datobs_all, datrnd_all] = cluster_test_helper(M_all_mahal, 10000,'t');
[h_all, p_all, ~] = cluster_test(datobs_all(time>=0,1),datrnd_all(time>=0,:), 0, 0.05,0.01);

figure
line('XData', [min(time) max(time)], 'YData', [0 0], 'LineStyle', '-', ...
    'LineWidth', 1, 'Color','k');
hold all
plot(time,mean(M_all_mahal,2),'c',time,mean(-M_long_mahal,2),'b')
scatter(time([false(1,length(time(time<0))) h_long']),ones(1,length(find(h_long==1))).*0.0005,'b','filled')
scatter(time([false(1,length(time(time<0))) h_all']),ones(1,length(find(h_all==1))).*0.0015,'c','filled')
ylim([-0.002 0.04]);
xlim([min(time) max(time)]);
pbaspect([1.1,1,1])
xlabel('Time relative to memory item')
ylabel('Multivaraite discrimination')
