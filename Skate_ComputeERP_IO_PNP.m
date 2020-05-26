
%CONDITIONS:
%preferred, clockwise - non-preffered, CCW

%%
ccc
%
exp = 'Skateboard';

% FULL LIST
%  subs = {'100' '101' '102'	'103'	'104'	'106'	'107'	'108'	'109' '110'...
%          '111' '112' '113'	'114'	'115'	'116'	'117'	'118'	'119' '120'...
%          '121' '122' '123'	'124'	'125'	'126'	'127'	'128'	'129' '130'...
%          '131' '132' '133'	'134' '135' '136' '137'};
% is_goofy = [0,	0,	0,	0,	0,	0,	1,	0,	1,	1,	0,	0,	1,	1,	1,	0,	0,	0,...
%             0,	0,	0,	1,	0,	0,	0,	0,	0,	1,	1,	0,	0,	0,	1,	1, 0, 1, 0];


% % %Full list minus subs that have loud FZ 
% subs = {'100'	'101'	'102'	'103'	'104'	'106'	'107'	'108'	'109'	'110'...
%     '111'	'112'	'113'	'114'	'115'	'116'	'117'	'119'	'120'...
%     '122'	'123'	'125'	'126'	'127'	'129'	'130'...
%     '131'	'132'	'133'	'134'	'135'	'136'	'137'};
% is_goofy = [0,	0,	0,	0,	0,	0,	1,	0,	1,	1,	0,	0,	1,	1,...
%     1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	0,	1,	1,    0,	1,	0];

subs = {'100'	'101'	'102'	'103'	'104'	'106'	'108'	'109'	'110'...
    '111'	'112'	'113'	'115'	'116'	'117'	'119'	'120'...
    '122'	'123'	'125'	'126'	'127'	'129'	...
    '131'	'132'	'133'	'134'	'136'	'137'};
is_goofy = [0,	0,	0,	0,	0,	0,	0,	1,	1,	0,	0,	1,	...
    1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	1,	1,    1,	0];

nsubs = length(subs); 
conds =  {'P_CW';'P_CCW'; 'NP_CW'; 'NP_CCW'};
conds_lab = {'Preferred Clockwise'; 'Preferred Counterclockwise'; 'Non-preferred Clockwise'; 'Non-preferred Counterclockwise'};
pref_lab = {'Preferred'; 'Non-Preferred'};
facing_lab = {'Facing Inside'; 'Facing Outside'};
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\Winter2019\'; %M:\Data\Skateboard\Winter2019
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%%
for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_Targets.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_JK\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = pop_loadset('filename',[Filename '_Standards.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_JK\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        
    end
end
eeglab redraw
%%
%subject erps
% electrode = 13;%this is PZ in this electrode map
%electrode = [13,14]; %to average for electrodes 

erp_out = [];
for i_sub = 1:nsubs
    fprintf(['Subject ' num2str(i_sub)])
    for i_cond = 1:nconds
        %average over trials (3rd dimension)
        erp_out(:,1,:,i_cond,i_sub) = mean(ALLEEG(1+ 2*((i_sub-1)*nconds+(i_cond-1))).data,3)'; %Targets
        erp_out(:,2,:,i_cond,i_sub) = mean(ALLEEG(2+ 2*((i_sub-1)*nconds+(i_cond-1))).data,3)'; %standards
    end
    if is_goofy(i_sub) == 1
        erp_out_face_in(:,:,:,i_sub) = squeeze(nanmean(erp_out(:,:,:,[2,3],i_sub),4));
        erp_out_face_out(:,:,:,i_sub) = squeeze(nanmean(erp_out(:,:,:,[1,4],i_sub),4));
    elseif is_goofy(i_sub) == 0
        erp_out_face_in(:,:,:,i_sub) = squeeze(nanmean(erp_out(:,:,:,[1,4],i_sub),4));
        erp_out_face_out(:,:,:,i_sub) = squeeze(nanmean(erp_out(:,:,:,[2,3],i_sub),4));
    end
     erp_out_p(:,:,:,i_sub) = squeeze(nanmean(erp_out(:,:,:,[1,2],i_sub),4)); %preferred
     erp_out_np(:,:,:,i_sub) = squeeze(nanmean(erp_out(:,:,:,[3,4],i_sub),4)); %non-preferred
     erp_out_CW(:,:,:,i_sub) = squeeze(nanmean(erp_out(:,:,:,[1,3],i_sub),4)); %clockwise merged
     erp_out_CCW(:,:,:,i_sub) = squeeze(nanmean(erp_out(:,:,:,[2,4],i_sub),4)); %counterclockwise merged 
end

%%
%grand average plots + difference
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
%time_window = find(EEG.times>-200,1)-1:find(EEG.times>1000,1)-2;
electrode = 15;
figure;
for i_cond = 1:nconds
    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'r';
        case 3
            colour = 'g';
        case 4
            colour = 'm';
    end
    
    subplot(2,nconds,i_cond);
    boundedline(EEG.times,squeeze(mean(erp_out(:,1,electrode,i_cond,:),5)),squeeze(std(erp_out(:,1,electrode,i_cond,:),[],5))./sqrt(nsubs),colour,...
        EEG.times,squeeze(mean(erp_out(:,2,electrode,i_cond,:),5)),squeeze(std(erp_out(:,2,electrode,i_cond,:),[],5))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-6 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(conds_lab{i_cond});
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    xlim([-200 1000]);
    legend({'Target', 'Standard'},'Location','northeast')
%     fill([195;195;295;295],[-6;12;12;-6],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
%     fill([375;375;525;525],[-6;12;12;-6],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
%     
    subplot(2,nconds,nconds+i_cond);
    boundedline(EEG.times,squeeze(mean(erp_diff_out(:,electrode,i_cond,:),4)),squeeze(std(erp_diff_out(:,electrode,i_cond,:),[],4))./sqrt(nsubs),colour);
   %boundedline(EEG.times,squeeze(mean(mean(erp_diff_out(:,electrode,i_cond,:),2),4)),squeeze(std(erp_diff_out(:,electrode,i_cond,:),[],4))./sqrt(nsubs),colour);
   ...for the averaging

    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-6 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Difference Wave');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    xlim([-200 1000]);
    legend({'Difference-Wave'},'Location','northeast');

%     fill([195;195;295;295],[-6;12;12;-6],'y','FaceAlpha',0.1, 'EdgeAlpha', '0.01');
%     fill([375;375;525;525],[-6;12;12;-6],'y','FaceAlpha',0.1,'EdgeAlpha', '0.01');
end
%%
%single plot for manuscript 

figure;
for i_cond = 1:nconds

    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'r';
        case 3
            colour = 'g';
        case 4
            colour = 'm';
    end
    
    electrode = 13;

    subplot(2,nconds,i_cond);
    hold on

    boundedline(EEG.times,squeeze(mean(erp_out(:,1,13,i_cond,:),5)),squeeze(std(erp_out(:,1,13,i_cond,:),[],5))./sqrt(nsubs),colour,...
        EEG.times,squeeze(mean(erp_out(:,2,13,i_cond,:),5)),squeeze(std(erp_out(:,2,13,i_cond,:),[],5))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-6 3]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
   % title(conds_lab{i_cond});
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
%     xlim([-200 1000]);
    legend({'Target', 'Standard'},'Location','northeast', 'Autoupdate', 'off');
    fill([195;195;295;295],[-6;3;3;-6],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
    %fill([375;375;525;525],[-6;2;2;-6],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
    hold off
   
    electrode = 15;
    subplot(2,nconds,nconds+i_cond);
    hold on
    boundedline(EEG.times,squeeze(mean(erp_out(:,1,electrode,i_cond,:),5)),squeeze(std(erp_out(:,1,electrode,i_cond,:),[],5))./sqrt(nsubs),colour,...
        EEG.times,squeeze(mean(erp_out(:,2,electrode,i_cond,:),5)),squeeze(std(erp_out(:,2,electrode,i_cond,:),[],5))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-4 6.5]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
   % title(conds_lab{i_cond});
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    xlim([-200 1000]);
    legend({'Target', 'Standard'},'Location','northeast','AutoUpdate','off');
    %fill([195;195;295;295],[-4;6.5;6.5;-4],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
    fill([350;350;525;525],[-4;6.5;6.5;-4],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
    hold off

end
%%
%Within subject anova - reran with JASP instead for simplicity 
% electrode = 15;
% time_window = find(EEG.times>195,1)-1:find(EEG.times>295,1)-2;
% time_window = find(EEG.times>350,1)-1:find(EEG.times>500,1)-2;
% 
% %p_targ = squeeze(mean(erp_out_p(time_window,1,electrode,:),1));
% p_cw_targ = squeeze(mean(erp_out(time_window,1,electrode,1,:),1));
% p_cw_stan = squeeze(mean(erp_out(time_window,2,electrode,1,:),1));
% p_ccw_targ = squeeze(mean(erp_out(time_window,1,electrode,2,:),1));
% p_ccw_stan = squeeze(mean(erp_out(time_window,2,electrode,2,:),1));
% np_cw_targ = squeeze(mean(erp_out(time_window,1,electrode,3,:),1));
% np_cw_stan = squeeze(mean(erp_out(time_window,2,electrode,3,:),1));
% np_ccw_targ = squeeze(mean(erp_out(time_window,1,electrode,4,:),1));
% np_ccw_stan = squeeze(mean(erp_out(time_window,2,electrode,4,:),1));
% 
% %p_targ_subj = 1:1:length(p_targ);
% p_cw_targ_subj = 1:1:length(p_cw_targ);
% p_cw_stan_subj = 1:1:length(p_cw_stan);
% p_ccw_targ_subj = 1:1:length(p_ccw_targ);
% p_ccw_stan_subj = 1:1:length(p_ccw_stan);
% np_cw_targ_subj = 1:1:length(np_cw_targ);
% np_cw_stan_subj = 1:1:length(np_cw_stan);
% np_ccw_targ_subj = 1:1:length(np_ccw_targ);
% np_ccw_stan_subj = 1:1:length(np_ccw_stan);
% 
% 
% %p_targ_tones = ones(length(p_targ),1)+0;
% p_cw_targ_tones = ones(length(p_cw_targ),1)+0;
% p_cw_stan_tones = ones(length(p_cw_stan),1)+1;
% p_ccw_targ_tones = ones(length(p_ccw_targ),1)+0;
% p_ccw_stan_tones = ones(length(p_ccw_stan),1)+1;
% np_cw_targ_tones = ones(length(np_cw_targ),1)+0;
% np_cw_stan_tones = ones(length(np_cw_stan),1)+1;
% np_ccw_targ_tones = ones(length(np_ccw_targ),1)+0;
% np_ccw_stan_tones = ones(length(np_ccw_stan),1)+1;
% 
% % p_targ_cond = ones(length(p_targ),1)+0;
% p_cw_targ_cond = ones(length(p_cw_targ),1)+0;
% p_cw_stan_cond = ones(length(p_cw_stan),1)+0;
% p_ccw_targ_cond = ones(length(p_ccw_targ),1)+1;
% p_ccw_stan_cond = ones(length(p_ccw_stan),1)+1;
% np_cw_targ_cond = ones(length(np_cw_targ),1)+2;
% np_cw_stan_cond = ones(length(np_cw_stan),1)+2;
% np_ccw_targ_cond = ones(length(np_ccw_targ),1)+3;
% np_ccw_stan_cond = ones(length(np_ccw_stan),1)+3;
% 
% 
% %       [p_targ,p_stan,np_targ,np_stan;...
% %     p_targ_subj,p_stan_subj,np_targ_subj,np_stan_subj;...
% %     p_targ_tones,p_stan_tones,np_targ_tones,np_stan_tones;...
% %     p_targ_cond,p_stan_cond,np_targ_cond,np_stan_cond]
% % [p_cw_targ,p_cw_stan,p_ccw_targ,p_ccw_stan,np_cw_targ,np_cw_stan,np_ccw_targ,np_ccw_stan;...
% % p_cw_targ_subj,p_cw_stan_subj,p_ccw_targ_subj,p_ccw_stan_subj,np_cw_targ_subj,np_cw_stan_subj,np_ccw_targ_subj,np_ccw_stan_subj;...
% % p_cw_targ_tones,p_cw_stan_tones,p_ccw_targ_tones,p_ccw_stan_tones,np_cw_targ_tones,np_cw_stan_tones,np_ccw_targ_tones,np_ccw_stan_tones;...
% % p_cw_targ_cond,p_cw_stan_cond,p_ccw_targ_cond,p_ccw_stan_cond,np_cw_targ_cond,np_cw_stan_cond,np_ccw_targ_cond,np_ccw_stan_cond];
% 
% % Y = [p_targ;p_stan;np_targ;np_stan];
% % S = [p_targ_subj';p_stan_subj';np_targ_subj';np_stan_subj'];
% % F1 = [p_targ_tones;p_stan_tones;np_targ_tones;np_stan_tones];
% % F2 = [p_targ_cond;p_stan_cond;np_targ_cond;np_stan_cond];
% % FACTNAMES = {'TONE TYPE','CONDITION'};
% % stats = rm_anova2(Y,S,F1,F2,FACTNAMES);
% Y = [p_cw_targ;p_cw_stan;p_ccw_targ;p_ccw_stan;np_cw_targ;np_cw_stan;np_ccw_targ;np_ccw_stan];
% S = [p_cw_targ_subj';p_cw_stan_subj';p_ccw_targ_subj';p_ccw_stan_subj';np_cw_targ_subj';np_cw_stan_subj';np_ccw_targ_subj';np_ccw_stan_subj'];
% F1 = [p_cw_targ_tones;p_cw_stan_tones;p_ccw_targ_tones;p_ccw_stan_tones;np_cw_targ_tones;np_cw_stan_tones;np_ccw_targ_tones;np_ccw_stan_tones];
% F2 = [p_cw_targ_cond;p_cw_stan_cond;p_ccw_targ_cond;p_ccw_stan_cond;np_cw_targ_cond;np_cw_stan_cond;np_ccw_targ_cond;np_ccw_stan_cond];
% FACTNAMES = {'TONE TYPE','CONDITION'};
% stats = rm_anova2(Y,S,F1,F2,FACTNAMES);
% 
% data_anova = [Y, S, F1, F2];
% stats2 = rmanova2(data_anova, 0.05,1,1);

%% 
% All conditions averaged to find early and late ERP windows 

%Peak to find time window 
electrode = 15;

figure
hold on
boundedline(EEG.times,squeeze(mean(mean(erp_diff_out(:,electrode,:,:),4),3)),squeeze(std(0))./sqrt(nsubs),'k');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; 
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Pz');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
 ylim ([-4 6]);
%  fill([195;195;295;295],[-4;6;6;-4],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
 fill([350;350;500;500],[-4;6;6;-4],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
hold off

time_window = find(EEG.times>195,1)-1:find(EEG.times>295,1)-2;
[h p ci stat] = ttest(squeeze(mean(mean(erp_diff_out(time_window,electrode,:,:),1),3)),0,.05,'both',1) 

%Topographies

%difference topographies across all conditions for peak plot: 
figure('Color',[1 1 1]);
    subplot(1,4,1);
    set(gca,'Color',[1 1 1]);
    temp = mean(mean(mean(erp_diff_out(time_window,:,:,:),3),4),1);
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
    %title(conds_lab{i_cond});
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
    
% 2x4 anova targ/standards, conds(4). 
%rm_anova2 () 
    
%%t-tests Original conditions
%ttest to compare Target-Standard conditions P3

electrode = 15;
time_window = find(EEG.times>375,1)-1:find(EEG.times>525,1)-2;
for i_cond = 1:nconds
[h p ci stat] = ttest(squeeze(mean(erp_out(time_window,1,electrode,i_cond,:),1)),squeeze(mean(erp_out(time_window,2,electrode,i_cond,:),1)),.01,'both',1) 
end 
electrode = 15;
time_window_p3 = find(EEG.times>375,1)-1:find(EEG.times>525,1)-2;
p3_dif = squeeze(mean(erp_out(time_window,1,electrode,4,:),1)) - squeeze(mean(erp_out(time_window,2,electrode,4,:),1));
[h p ci stat] = ttest(p3_dif,0,.01,'left',1) 
mean (p3_dif)

%ttest to compare Target-Standard conditions MMN at Pz 
electrode = 15;
time_window_mmn = find(EEG.times>195,1)-1:find(EEG.times>295,1)-2;
mmn_dif = squeeze(mean(erp_out(time_window,1,electrode,4,:),1)) - squeeze(mean(erp_out(time_window,2,electrode,4,:),1));
[h p ci stat] = ttest(mmn_dif,0,.05,'left',1) 
mean (mmn_dif)

%ttest to compare Target-Standard conditions MMN
electrode = 13;
time_window = find(EEG.times>195,1)-1:find(EEG.times>295,1)-2;
mmn_dif = squeeze(mean(erp_out(time_window,1,electrode,4,:),1)) - squeeze(mean(erp_out(time_window,2,electrode,4,:),1));
[h p ci stat] = ttest(mmn_dif,0,.05,'both',1) 
mean (mmn_dif)

%t-test original conds p3
electrode = 15;
time_window = find(EEG.times>350,1)-1:find(EEG.times>500,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,electrode,3,:),1)),squeeze(mean(erp_diff_out(time_window,electrode,4,:),1)),.05,'both',1) 

%t-test original conds MMN
electrode = 13;
time_window = find(EEG.times>200,1)-1:find(EEG.times>300,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,electrode,2,:),1)),squeeze(mean(erp_diff_out(time_window,electrode,4,:),1)),.05,'both',1) 


%Topographies

%difference topographies
%time_window = find(EEG.times>200,1)-1:find(EEG.times>300,1)-2;
figure('Color',[1 1 1]);
for i_cond = 1:nconds
    subplot(1,nconds,i_cond);
    set(gca,'Color',[1 1 1]);
    temp = mean(mean(erp_diff_out(time_window_mmn,:,i_cond,:),4),1);
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
    title(conds_lab{i_cond});
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
    
end

%%
%                   by preference
%
 %pref-nonpref + difference waves PZ and FZ 
 
erp_diff_out_p = squeeze(erp_out_p(:,1,:,:)-erp_out_p(:,2,:,:));
erp_diff_out_np = squeeze(erp_out_np(:,1,:,:)-erp_out_np(:,2,:,:));

electrode = 13;
clear L
figure;
hold on
subplot(2,4,1);
    boundedline(EEG.times,squeeze(mean(erp_out_p(:,1,electrode,:),4)),squeeze(std(erp_out_p(:,1,electrode,:),[],4))./sqrt(nsubs),'r',...
        EEG.times,squeeze(mean(erp_out_p(:,2,electrode,:),4)),squeeze(std(erp_out_p(:,2,electrode,:),[],4))./sqrt(nsubs),'k');   
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-6 3]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    %title('Preferred');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    legend({'Target', 'Standard'},'Location','northeast')
    xlim ([-200 1000]);
    
    subplot(2,4,2);
    boundedline(EEG.times,squeeze(mean(erp_out_np(:,1,electrode,:),4)),squeeze(std(erp_out_np(:,1,electrode,:),[],4))./sqrt(nsubs),'b',...
        EEG.times,squeeze(mean(erp_out_np(:,2,electrode,:),4)),squeeze(std(erp_out_np(:,2,electrode,:),[],4))./sqrt(nsubs),'k');    
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-6 3]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    %title('Non-preferred');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)'); 
    legend({'Target', 'Standard'},'Location','northeast')
    xlim ([-200 1000]);

    
    subplot(2,4,3);
    boundedline(EEG.times,squeeze(mean(erp_diff_out_p(:,electrode,:),3)),squeeze(std(erp_diff_out_p(:,electrode,:),[],3))./sqrt(nsubs),'b',... %erp_diff_out_p(:,electrode,:),[],3
        EEG.times,squeeze(mean(erp_diff_out_np(:,electrode,:),3)),squeeze(std(erp_diff_out_np(:,electrode,:),[],3))./sqrt(nsubs),'r'); %erp_diff_out_np(:,2,electrode,i_cond,:),[],4
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-6 3]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    %title('Difference-Wave by Preference');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    xlim ([-200 1000]);
    fill([195;195;295;295],[-6;3;3;-6],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
    L(1) = plot(nan, nan, 'b');
    L(2) = plot(nan, nan, 'r');
    legend(L, {'Difference-wave preferred', 'Difference-wave non-preferred'},'Location','northeast')
    hold off
    
    electrode = 15;
    subplot(2,4,5);
    boundedline(EEG.times,squeeze(mean(erp_out_p(:,1,electrode,:),4)),squeeze(std(erp_out_p(:,1,electrode,:),[],4))./sqrt(nsubs),'r',...
        EEG.times,squeeze(mean(erp_out_p(:,2,electrode,:),4)),squeeze(std(erp_out_p(:,2,electrode,:),[],4))./sqrt(nsubs),'k');   
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-4 6.5]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    %title('Preferred');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    legend({'Target', 'Standard'},'Location','northeast')
    xlim ([-200 1000]);
    
    subplot(2,4,6);
    boundedline(EEG.times,squeeze(mean(erp_out_np(:,1,electrode,:),4)),squeeze(std(erp_out_np(:,1,electrode,:),[],4))./sqrt(nsubs),'b',...
        EEG.times,squeeze(mean(erp_out_np(:,2,electrode,:),4)),squeeze(std(erp_out_np(:,2,electrode,:),[],4))./sqrt(nsubs),'k');    
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-4 6.5]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    %title('Non-preferred');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)'); 
    legend({'Target', 'Standard'},'Location','northeast')
    xlim ([-200 1000]);

    clear L
    subplot(2,4,7);
    boundedline(EEG.times,squeeze(mean(erp_diff_out_p(:,electrode,:),3)),squeeze(std(erp_diff_out_p(:,electrode,:),[],3))./sqrt(nsubs),'b',... %erp_diff_out_p(:,electrode,:),[],3
        EEG.times,squeeze(mean(erp_diff_out_np(:,electrode,:),3)),squeeze(std(erp_diff_out_np(:,electrode,:),[],3))./sqrt(nsubs),'r'); %erp_diff_out_np(:,2,electrode,i_cond,:),[],4
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-4 6.5]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    %title('Difference-Wave by Preference');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
%     fill([195;195;295;295],[-4;6.5;6.5;-4],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
    fill([375;375;525;525],[-4;6.5;6.5;-4],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
    L(1) = plot(nan, nan, 'b');
    L(2) = plot(nan, nan, 'r');
    legend(L, {'Difference-wave preferred', 'Difference-wave non-preferred'},'Location','northeast')
    hold off
%     

figure('Color',[1 1 1]);
% for i_cond = 1:nconds
    subplot(1,2,1);
    set(gca,'Color',[1 1 1]);
    temp = mean(mean(erp_diff_out_p(time_window_p3,:,:),3),1);
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
    title('Preferred');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
    
    subplot(1,2,2);
    set(gca,'Color',[1 1 1]);
    temp = mean(mean(erp_diff_out_np(time_window_p3,:,:),3),1);
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
    title('Non-Preferred');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
    

    
    
%%
    %ttests by preference
    
erp_diff_out_p = squeeze(erp_out_p(:,1,:,:)-erp_out_p(:,2,:,:));

erp_diff_out_np = squeeze(erp_out_np(:,1,:,:)-erp_out_np(:,2,:,:));

%t-test original conds p3 
electrode = 15;
time_window = find(EEG.times>350,1)-1:find(EEG.times>500,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_diff_out_p(time_window,electrode,:),1)),squeeze(mean(erp_diff_out_np(time_window,electrode,:),1)),.05,'both',1) 
p3_dif_ttest = squeeze(mean(erp_diff_out_p(time_window,electrode,:),1)) - squeeze(mean(erp_diff_out_np(time_window,electrode,:),1))
[h p ci stat] = ttest(p3_dif_ttest,0,.05,'right',1) 
mean (p3_dif_ttest)

%t-test preference MMN
electrode = 13;
time_window = find(EEG.times>195,1)-1:find(EEG.times>295,1)-2;
mmn_dif_pnp = squeeze(mean(erp_diff_out_p(time_window,electrode,:),1)) - squeeze(mean(erp_diff_out_np(time_window,electrode,:),1));
[h p ci stat] = ttest(squeeze(mean(erp_diff_out_p(time_window,electrode,:),1)),squeeze(mean(erp_diff_out_np(time_window,electrode,:),1)),.05,'both',1) 
[h p ci stat] = ttest(mmn_dif_pnp,0,.05,'both',1)
mean (mmn_dif_pnp)


%Repeated measures anova Preference/orientation 
electrode = 13;
time_window = find(EEG.times>195,1)-1:find(EEG.times>295,1)-2;
p_targ = squeeze(mean(erp_out_p(time_window,1,electrode,:),1));   
p_stan = squeeze(mean(erp_out_p(time_window,2,electrode,:),1));
np_targ = squeeze(mean(erp_out_np(time_window,1,electrode,:),1));   
np_stan = squeeze(mean(erp_out_np(time_window,2,electrode,:),1));

p_targ_subj = 1:1:length(p_targ);
p_stan_subj = 1:1:length(p_stan);
np_targ_subj = 1:1:length(np_targ);
np_stan_subj = 1:1:length(np_stan);

p_targ_tones = ones(length(p_targ),1)+0;
p_stan_tones = ones(length(p_stan),1)+1;
np_targ_tones = ones(length(np_targ),1)+0;
np_stan_tones = ones(length(np_stan),1)+1;

p_targ_cond = ones(length(p_targ),1)+0;
p_stan_cond = ones(length(p_stan),1)+0;
np_targ_cond = ones(length(np_targ),1)+1;
np_stan_cond = ones(length(np_stan),1)+1;

% [p_targ,p_stan,np_targ,np_stan;...
%     p_targ_subj,p_stan_subj,np_targ_subj,np_stan_subj;...
%     p_targ_tones,p_stan_tones,np_targ_tones,np_stan_tones;...
%     p_targ_cond,p_stan_cond,np_targ_cond,np_stan_cond]

Y = [p_targ;p_stan;np_targ;np_stan];
S = [p_targ_subj';p_stan_subj';np_targ_subj';np_stan_subj'];
F1 = [p_targ_tones;p_stan_tones;np_targ_tones;np_stan_tones];
F2 = [p_targ_cond;p_stan_cond;np_targ_cond;np_stan_cond];
FACTNAMES = {'TONE TYPE','CONDITION'};
statspref = rm_anova2(Y,S,F1,F2,FACTNAMES);


data_anova = [Y, S, F1, F2];
statspref2 = rmanova2(data_anova, 0.05,1,1);

% ranovatbl = ranova (

% TOPOGRAPHIES BY PREFERENCE 
%time_window = find(EEG.times>375,1)-1:find(EEG.times>525,1)-2;
time_window = find(EEG.times>375,1)-1:find(EEG.times>575,1)-2;

figure('Color',[1 1 1]);
subplot(1,2,1)
    set(gca,'Color',[1 1 1]);
    temp = mean(erp_diff_out_p(time_window_p3,:,:),1);
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
    title('Preferred Stance');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
    
    subplot(1,2,2)  
    temp = mean(erp_diff_out_np(time_window_p3,:,:),1);
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
    title('Non-Preferred Stance');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
%%    
%Targets on same axis by preference
figure;
boundedline(EEG.times,squeeze(mean(erp_out_p(:,1,electrode,:),4)),squeeze(std(erp_out_p(:,1,electrode,:),[],4))./sqrt(nsubs),'r',...
        EEG.times,squeeze(mean(erp_out_np(:,1,electrode,:),4)),squeeze(std(erp_out_np(:,1,electrode,:),[],4))./sqrt(nsubs),'b');  
    
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Target Tones');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
L(1) = plot(nan, nan, 'r');
L(2) = plot(nan, nan, 'b');
legend(L, {'preferred stance', 'non-preferred stance'},'Location','northeast', 'Autoupdate', 'off')
fill([125;125;225;225],[-10;12;12;-10],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
fill([225;225;325;325],[-10;12;12;-10],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');

%ttest for targets by preference at N1
time_window = find(EEG.times>125,1)-1:find(EEG.times>200,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_p(time_window,1,electrode,:),1)),squeeze(mean(erp_out_np(time_window,1,electrode,:),1)),.05,'both',1) 

%ttest for targets by preference at P2
time_window = find(EEG.times>225,1)-1:find(EEG.times>325,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_p(time_window,1,electrode,:),1)),squeeze(mean(erp_out_np(time_window,1,electrode,:),1)),.05,'both',1) 

%Standards on same axis by preference
figure;
boundedline(EEG.times,squeeze(mean(erp_out_p(:,2,electrode,:),4)),squeeze(std(erp_out_p(:,2,electrode,:),[],4))./sqrt(nsubs),'r',...
        EEG.times,squeeze(mean(erp_out_np(:,2,electrode,:),4)),squeeze(std(erp_out_np(:,2,electrode,:),[],4))./sqrt(nsubs),'b');  
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-10 15],'color','k');
title('Standard Tones');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
fill([125;125;225;225],[-10;12;12;-10],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
fill([225;225;325;325],[-10;12;12;-10],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');

L(1) = plot(nan, nan, 'r');
L(2) = plot(nan, nan, 'b');
legend(L, {'preferred stance', 'non-preferred stance'},'Location','northeast')


%ttest for targets by preference at N1
time_window = find(EEG.times>125,1)-1:find(EEG.times>225,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_p(time_window,2,electrode,:),1)),squeeze(mean(erp_out_np(time_window,2,electrode,:),1)),.05,'both',1) 

%ttest for targets by preference at P2
time_window = find(EEG.times>225,1)-1:find(EEG.times>325,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_p(time_window,2,electrode,:),1)),squeeze(mean(erp_out_np(time_window,2,electrode,:),1)),.05,'both',1) 


%%
%fOUR CONDS TARG STANDS SAME AXIS 
%Targets on same axis by preference
electrode = 15;
figure;
for i_cond = 1:nconds
    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'r';
        case 3
            colour = 'g';
        case 4
            colour = 'm';
    end
    
%     subplot(2,nconds,i_cond);
    boundedline(EEG.times,squeeze(mean(erp_out(:,1,electrode,i_cond,:),5)),squeeze(std(0))./sqrt(nsubs),colour);,...
%         EEG.times,squeeze(mean(erp_out(:,2,electrode,i_cond,:),5)),squeeze(std(erp_out(:,2,electrode,i_cond,:),[],5))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-6 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Target Tones');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    xlim([-200 1000]);
    
fill([125;125;200;200],[-10;12;12;-10],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
fill([200;200;300;300],[-10;12;12;-10],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
   
end

L(1) = plot(nan, nan, 'b');
L(2) = plot(nan, nan, 'r');
L(3) = plot(nan, nan, 'g');
L(4) = plot(nan, nan, 'm');
legend(L, {'preferred Clockwise','preferred Counterclockwise', 'non-preferred Clockwise','Non-preferred Clockwise', 'Non-preferred Counterclockwise'},'Location','northeast')

figure;
for i_cond = 1:nconds
    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'r';
        case 3
            colour = 'g';
        case 4
            colour = 'm';
    end
    
%     subplot(2,nconds,i_cond);
    boundedline(EEG.times,squeeze(mean(erp_out(:,2,electrode,i_cond,:),5)),squeeze(std(0))./sqrt(nsubs),colour);
%         EEG.times,squeeze(mean(erp_out(:,2,electrode,i_cond,:),5)),squeeze(std(erp_out(:,2,electrode,i_cond,:),[],5))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-6 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Standard Tones');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    xlim([-200 1000]);
    
fill([125;125;200;200],[-10;12;12;-10],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
fill([200;200;300;300],[-10;12;12;-10],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
    
end

L(1) = plot(nan, nan, 'b');
L(2) = plot(nan, nan, 'r');
L(3) = plot(nan, nan, 'g');
L(4) = plot(nan, nan, 'm');
legend(L, {'preferred Clockwise','preferred Counterclockwise', 'non-preferred Clockwise','Non-preferred Clockwise', 'Non-preferred Counterclockwise'},'Location','northeast')

%%
%peaks n1 p2
%Peak to find time window 
electrode = 15;

figure
hold on
    boundedline(EEG.times,squeeze(mean(mean(erp_out(:,2,electrode,i_cond,:),4),5)),squeeze(std(0))./sqrt(nsubs),colour);
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; 
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title(electrode);
xlabel('Time (ms)');
ylabel('Voltage (uV)');
 ylim ([-6 6]);
%  fill([195;195;295;295],[-4;6;6;-4],'w','FaceAlpha',0.1, 'EdgeAlpha', '1', 'Linestyle', ':');
 %fill([350;350;500;500],[-4;6;6;-4],'w','FaceAlpha',0.1,'EdgeAlpha', '1', 'Linestyle', ':');
hold off

%ttest for targets by preference at N1
time_window = find(EEG.times>125,1)-1:find(EEG.times>200,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out(time_window,1,electrode,1,:),1)),squeeze(mean(erp_out_np(time_window,1,electrode,:),1)),.05,'both',1) 

%ttest for targets by preference at P2
time_window = find(EEG.times>200,1)-1:find(EEG.times>300,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_np(time_window,1,electrode,:),1)),squeeze(mean(erp_out_np(time_window,1,electrode,:),1)),.05,'both',1) 

%%
    %Peak to find time window 
mean_facein = squeeze(mean(erp_diff_out_face_in(:,electrode,:),3));
mean_faceout = squeeze(mean(erp_diff_out_face_out(:,electrode,:),3));
average_IO = mean_facein+mean_faceout/2;
figure
boundedline(EEG.times,average_PNP,squeeze(std(0))./sqrt(nsubs),colour);
 
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-6 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Difference Wave');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    xlim ([-200 1000]);
    ylim ([-8 12]);
    %%
%difference topographies by preference
%time_window = find(EEG.times>375,1)-1:find(EEG.times>525,1)-2;
time_window = find(EEG.times>195,1)-1:find(EEG.times>295,1)-2;

figure('Color',[1 1 1]);
subplot(1,2,1)
    set(gca,'Color',[1 1 1]);
    temp = mean(erp_diff_out_face_in(time_window,:,:),1);
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
    title('Face In');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
    
    subplot(1,2,2)  
    temp = mean(erp_diff_out_face_out(time_window,:,:),1);
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4])
    title('Face Out');
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');

%%
%difference on same axis original conditions
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
figure   
    boundedline(EEG.times,squeeze(mean(erp_diff_out(:,electrode,1,:),4)),squeeze(std(erp_diff_out(:,electrode,1,:),[],4))./sqrt(nsubs),'b',...
        EEG.times,squeeze(mean(erp_diff_out(:,electrode,2,:),4)),squeeze(std(erp_diff_out(:,electrode,2,:),[],4))./sqrt(nsubs),'r',...
        EEG.times,squeeze(mean(erp_diff_out(:,electrode,3,:),4)),squeeze(std(erp_diff_out(:,electrode,3,:),[],4))./sqrt(nsubs),'g',...
        EEG.times,squeeze(mean(erp_diff_out(:,electrode,4,:),4)),squeeze(std(erp_diff_out(:,electrode,4,:),[],4))./sqrt(nsubs),'m');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Difference-Wave');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r');
L(3) = plot(nan, nan, 'g');
L(4) = plot(nan, nan, 'm');
legend(L, {'Preferred Clockwise'; 'Preferred Counterclockwise'; 'Non-preferred Clockwise'; 'Non-preferred Counterclockwise'},'Location','northwest')

%comparing non-preferred (original) conditions clockwise vs counterclock wise 
figure   
    boundedline(EEG.times,squeeze(mean(erp_diff_out(:,electrode,3,:),4)),squeeze(std(erp_diff_out(:,electrode,3,:),[],4))./sqrt(nsubs),'g',...
        EEG.times,squeeze(mean(erp_diff_out(:,electrode,4,:),4)),squeeze(std(erp_diff_out(:,electrode,4,:),[],4))./sqrt(nsubs),'m');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Difference-Wave');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    
L(1) = plot(nan, nan, 'g');
L(2) = plot(nan, nan, 'm');
legend(L, {'Non-preferred Clockwise'; 'Non-preferred Counterclockwise'},'Location','northwest');


%ERPs by preference/orientation same axis
figure
    boundedline(EEG.times,squeeze(mean(erp_diff_out_p(:,electrode,:),3)),squeeze(std(erp_diff_out_p(:,electrode,:),[],3))./sqrt(nsubs),'b',... %erp_diff_out_p(:,electrode,:),[],3
        EEG.times,squeeze(mean(erp_diff_out_np(:,electrode,:),3)),squeeze(std(erp_diff_out_np(:,electrode,:),[],3))./sqrt(nsubs),'r',... %erp_diff_out_np(:,2,electrode,i_cond,:),[],4
        EEG.times,squeeze(mean(erp_diff_out_face_in(:,electrode,:),3)),squeeze(std(erp_diff_out_face_in(:,electrode,:),[],3))./sqrt(nsubs),'g',... %erp_diff_out_p(:,electrode,:),[],3
        EEG.times,squeeze(mean(erp_diff_out_face_out(:,electrode,:),3)),squeeze(std(erp_diff_out_face_out(:,electrode,:),[],3))./sqrt(nsubs),'m'); %erp_diff_out_np(:,2,electrode,i_cond,:),[],4   
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Difference-Waves');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r');
L(3) = plot(nan, nan, 'g');
L(4) = plot(nan, nan, 'm');
legend(L, {'preferred stance', 'non-preferred stance','Facing Inside', 'Facing Outside'},'Location','northwest')
    
%%
%CLOCKWISE VS COUNTERCLOCKWISE orientation


figure
subplot (1,3,1)
boundedline(EEG.times,squeeze(mean(erp_out_CW(:,2,electrode,:),4)),squeeze(std(erp_out_CW(:,2,electrode,:),[],4))./sqrt(nsubs),'k',...
    EEG.times,squeeze(mean(erp_out_CW(:,1,electrode,:),4)),squeeze(std(erp_out_CW(:,1,electrode,:),[],4))./sqrt(nsubs),'b');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Clockwise');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
legend('Standards', 'Targets','Location','northwest')
xlim([-200 1000])

subplot (1,3,2)
boundedline(EEG.times,squeeze(mean(erp_out_CCW(:,2,electrode,:),4)),squeeze(std(erp_out_CCW(:,2,electrode,:),[],4))./sqrt(nsubs),'k',...
    EEG.times,squeeze(mean(erp_out_CCW(:,1,electrode,:),4)),squeeze(std(erp_out_CCW(:,1,electrode,:),[],4))./sqrt(nsubs),'g');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Counterclock-wise');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
legend('Standards', 'Targets','Location','northwest')
xlim([-200 1000])

subplot (1,3,3)
boundedline(EEG.times,squeeze(mean(erp_diff_out_CW(:,electrode,:),3)),squeeze(std(erp_diff_out_CW(:,electrode,:),[],3))./sqrt(nsubs),'r',...
    EEG.times,squeeze(mean(erp_diff_out_CCW(:,electrode,:),3)),squeeze(std(erp_diff_out_CCW(:,electrode,:),[],3))./sqrt(nsubs),'b');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('CW CCW Difference-Wave');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
legend('Clockwise', 'Counterclock-wise','Location','northwest')
xlim([-200 1000])
%%
%CLOCKWISE VS COUNTERCLOCKWISE STANDARDS AND TARGETS ON SAME AXIS


figure
subplot (1,2,1)
boundedline(EEG.times,squeeze(mean(erp_out_CW(:,1,electrode,:),4)),squeeze(std(erp_out_CW(:,1,electrode,:),[],4))./sqrt(nsubs),'k',...
    EEG.times,squeeze(mean(erp_out_CCW(:,1,electrode,:),4)),squeeze(std(erp_out_CCW(:,1,electrode,:),[],4))./sqrt(nsubs),'b');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Target Tone');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
legend('Clockwise', 'Counterclock-wise','Location','northwest')

subplot (1,2,2)
boundedline(EEG.times,squeeze(mean(erp_out_CW(:,2,electrode,:),4)),squeeze(std(erp_out_CW(:,2,electrode,:),[],4))./sqrt(nsubs),'r',...
    EEG.times,squeeze(mean(erp_out_CCW(:,2,electrode,:),4)),squeeze(std(erp_out_CCW(:,2,electrode,:),[],4))./sqrt(nsubs),'g');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Standard Tone');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
legend('Clockwise', 'Counterclock-wise','Location','northwest')

%t-test for targets N1 by clock orientation
time_window = find(EEG.times>125,1)-1:find(EEG.times>225,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_CW(time_window,1,electrode,:),1)),squeeze(mean(erp_out_CCW(time_window,1,electrode,:),1)),.05,'both',1) 

%ttest for targets by preference at P2
time_window = find(EEG.times>200,1)-1:find(EEG.times>300,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_CW(time_window,1,electrode,:),1)),squeeze(mean(erp_out_CCW(time_window,1,electrode,:),1)),.05,'both',1) 

%comparing clock orientation in difference waves 
erp_diff_out_CW = squeeze(erp_out_CW(:,1,:,:)-erp_out_CW(:,2,:,:));
erp_diff_out_CCW = squeeze(erp_out_CCW(:,1,:,:)-erp_out_CCW(:,2,:,:));

figure;
boundedline(EEG.times,squeeze(mean(erp_diff_out_CW(:,electrode,:),3)),squeeze(std(erp_diff_out_CW(:,electrode,:),[],3))./sqrt(nsubs),'r',...
    EEG.times,squeeze(mean(erp_diff_out_CCW(:,electrode,:),3)),squeeze(std(erp_diff_out_CCW(:,electrode,:),[],3))./sqrt(nsubs),'b');
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Clockwise vs CCW at Fz');
xlabel('Time (ms)');
ylabel('Voltage (uV)');
legend('Clockwise', 'Counterclock-wise','Location','northwest')

% TESTING FOR DIFFERENCES
time_window = find(EEG.times>300,1)-1:find(EEG.times>500,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_diff_out_CW(time_window,electrode,:),1)),squeeze(mean(erp_diff_out_CCW(time_window,electrode,:),1)),.05,'both',1) 

%find peak
figure;
boundedline(EEG.times,squeeze(mean(mean(erp_diff_out(:,electrode,:,:),3),4)), squeeze(std(0))./sqrt(nsubs),'k'),...
    
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-2.5 8],'color','k');
title('Target Grand ERPs');
xlabel('Time (ms)');
ylabel('Voltage (uV)');

%%
% TESTING FOR DIFFERENCES
%erp_diff_out(:,electrode,i_cond,:)

time_window = find(EEG.times>375,1)-1:find(EEG.times>575,1)-2;
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));

%%%timepoints X events X electrodes X conditions X participants%%%%

%  COMPARING DIFFERENCE WAVES for non-preferred conditions
[h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,electrode,3,:),1)),squeeze(mean(erp_diff_out(time_window,electrode,4,:),1)),.05,'both',1) 
%%

%SINGLE ERPS FOR INDIVIDUAL COPYING

erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
figure('Color',[1 1 1]);
for i_cond = 1:nconds
    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'r';
        case 3
            colour = 'g';
        case 4
            colour = 'm';
    end
    
    figure;
    boundedline(EEG.times,squeeze(mean(erp_out(:,1,electrode,i_cond,:),5)),squeeze(std(erp_out(:,1,electrode,i_cond,:),[],5))./sqrt(nsubs),colour,...
        EEG.times,squeeze(mean(erp_out(:,2,electrode,i_cond,:),5)),squeeze(std(erp_out(:,2,electrode,i_cond,:),[],5))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    if i_cond == 2
        legend('Targets','Standards','Location','NorthEast');
    elseif i_cond == 1
        legend('Targets','Standards','Location','NorthEast');
    end
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(conds_lab{i_cond});
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
end


