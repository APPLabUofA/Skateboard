
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


% %Full list minus subs that have loud FZ 
subs = {'100'	'101'	'102'	'103'	'104'	'106'	'107'	'108'	'109'	'110'...
    '111'	'112'	'113'	'114'	'115'	'116'	'117'	'119'	'120'...
    '122'	'123'	'125'	'126'	'127'	'129'	'130'...
    '131'	'132'	'133'	'134'	'135'	'136'	'137'};
is_goofy = [0,	0,	0,	0,	0,	0,	1,	0,	1,	1,	0,	0,	1,	1,...
    1,	0,	0,	0,	0,	1,	0,	0,	0,	0,	1,	0,	0,	0,	1,	1,    0,	1,	0];

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
        EEG = pop_loadset('filename',[Filename '_Corrected_Target.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_test\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = pop_loadset('filename',[Filename '_Corrected_Standard.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_test\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        
    end
end
eeglab redraw
%%
%subject erps
electrode = 15;%this is PZ in this electrode map
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
    
end

%individual subs
%  
% for i_sub = 1:nsubs
%     
%          figure
%         boundedline(EEG.times,squeeze(erp_out(:,1,electrode,:,i_sub)),squeeze(std(0))./sqrt(nsubs),'r',...
%             EEG.times,squeeze(erp_out(:,2,electrode,:,i_sub)),squeeze(std(0))./sqrt(nsubs),'k');
%         set(gca,'Color',[1 1 1]);
%         set(gca,'YDir','reverse');
%         axis tight; ylim([-20 20]);
%         line([-200 1000],[0 0],'color','k');
%         line([0 0],[-2.5 8],'color','k');
%         title(subs{i_sub});
%         xlabel('Time (ms)');
%         ylabel('Voltage (uV)');
%  end


%%
%                   by preference
%
 %pref-nonpref + difference waves
erp_diff_out_p = squeeze(erp_out_p(:,1,:,:)-erp_out_p(:,2,:,:));
erp_diff_out_np = squeeze(erp_out_np(:,1,:,:)-erp_out_np(:,2,:,:));

figure;
subplot(1,3,1);
    boundedline(EEG.times,squeeze(mean(erp_out_p(:,1,electrode,:),4)),squeeze(std(erp_out_p(:,1,electrode,:),[],4))./sqrt(nsubs),'r',...
        EEG.times,squeeze(mean(erp_out_p(:,2,electrode,:),4)),squeeze(std(erp_out_p(:,2,electrode,:),[],4))./sqrt(nsubs),'k');   
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Preferred');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    legend({'Target', 'Standard'},'Location','northwest')

subplot(1,3,2);
    boundedline(EEG.times,squeeze(mean(erp_out_np(:,1,electrode,:),4)),squeeze(std(erp_out_np(:,1,electrode,:),[],4))./sqrt(nsubs),'b',...
        EEG.times,squeeze(mean(erp_out_np(:,2,electrode,:),4)),squeeze(std(erp_out_np(:,2,electrode,i_cond,:),[],4))./sqrt(nsubs),'k');    
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Non-preferred');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)'); 
    legend({'Target', 'Standard'},'Location','northwest')
   
subplot(1,3,3);
    boundedline(EEG.times,squeeze(mean(erp_diff_out_p(:,electrode,:),3)),squeeze(std(erp_diff_out_p(:,electrode,:),[],3))./sqrt(nsubs),'b',... %erp_diff_out_p(:,electrode,:),[],3
        EEG.times,squeeze(mean(erp_diff_out_np(:,electrode,:),3)),squeeze(std(erp_diff_out_np(:,electrode,:),[],3))./sqrt(nsubs),'r'); %erp_diff_out_np(:,2,electrode,i_cond,:),[],4
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Difference-Wave by Preference');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    legend({'preferred stance', 'non-preferred stance'},'Location','northwest')

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
legend(L, {'preferred stance', 'non-preferred stance'},'Location','northeast')

%ttest for targets by preference at N1
time_window = find(EEG.times>125,1)-1:find(EEG.times>225,1)-2;
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
%  %Averaging between 2 electrodes
%  electrode = [13, 14];
% 
%  figure;
%     boundedline(EEG.times,squeeze(mean(mean(erp_out_p(:,1,electrode,:),3),4)),squeeze(std(erp_out_p(:,1,electrode,:),[],4))./sqrt(nsubs),'r',...
%         EEG.times,squeeze(mean(mean(erp_out_p(:,2,electrode,:),3),4)),squeeze(std(erp_out_p(:,2,electrode,:),[],4))./sqrt(nsubs),'k');    set(gca,'Color',[1 1 1]);
%     set(gca,'YDir','reverse');
%     axis tight; ylim([-8 12]);
%     line([-200 1000],[0 0],'color','k');
%     line([0 0],[-2.5 8],'color','k');
%     title('Preferred');
%     xlabel('Time (ms)');
%     ylabel('Voltage (uV)');
%     legend({'Target', 'Standard'},'Location','northwest')
% 
% figure;
% boundedline(EEG.times,squeeze(mean(mean(erp_out_np(:,1,electrode,:),3),4)),squeeze(std(erp_out_np(:,1,electrode,:),[],4))./sqrt(nsubs),'b',...
%         EEG.times,squeeze(mean(mean(erp_out_p(:,2,electrode,:),3),4)),squeeze(std(erp_out_p(:,2,electrode,i_cond,:),[],4))./sqrt(nsubs),'k');    
%     set(gca,'Color',[1 1 1]);
%     set(gca,'YDir','reverse');
%     axis tight; ylim([-8 12]);
%     line([-200 1000],[0 0],'color','k');
%     line([0 0],[-2.5 8],'color','k');
%     title('Non-preferred');
%     xlabel('Time (ms)');
%     ylabel('Voltage (uV)'); 
%     legend({'Target', 'Standard'},'Location','northwest')
%    
%     erp_diff_out_p = squeeze(erp_out_p(:,1,:,:)-erp_out_p(:,2,:,:));
% erp_diff_out_np = squeeze(erp_out_np(:,1,:,:)-erp_out_np(:,2,:,:));
% 
% figure;
% boundedline(EEG.times,squeeze(mean(mean(erp_diff_out_p(:,electrode,:),2),3)),squeeze(std(0))./sqrt(nsubs),'b',... %erp_diff_out_p(:,electrode,:),[],3
%         EEG.times,squeeze(mean(mean(erp_diff_out_np(:,electrode,:),2),3)),squeeze(std(0))./sqrt(nsubs),'r'); %erp_diff_out_np(:,electrode,:),[],3
%     set(gca,'Color',[1 1 1]);
%     set(gca,'YDir','reverse');
%     axis tight; ylim([-8 12]);
%     line([-200 1000],[0 0],'color','k');
%     line([0 0],[-2.5 8],'color','k');
%     title('Difference-Wave by Preference');
%     xlabel('Time (ms)');
%     ylabel('Voltage (uV)');
%     legend({'preferred stance', 'non-preferred stance'},'Location','northwest')

    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%
%face inside vs facing outside
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
subplot(1,3,1);
    boundedline(EEG.times,squeeze(mean(erp_out_face_in(:,1,electrode,:),4)),squeeze(std(erp_out_face_in(:,1,electrode,:),[],4))./sqrt(nsubs),'m',...
        EEG.times,squeeze(mean(erp_out_face_in(:,2,electrode,:),4)),squeeze(std(erp_out_face_in(:,2,electrode,:),[],4))./sqrt(nsubs),'k');    
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Facing in');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
     legend({'Target', 'Standard'},'Location','northwest')

 
subplot(1,3,2);

    boundedline(EEG.times,squeeze(mean(erp_out_face_out(:,1,electrode,:),4)),squeeze(std(erp_out_face_out(:,1,electrode,:),[],4))./sqrt(nsubs),'b',...
        EEG.times,squeeze(mean(erp_out_face_out(:,2,electrode,:),4)),squeeze(std(erp_out_face_out(:,2,electrode,:),[],4))./sqrt(nsubs),'k');    
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Facing out');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)'); 
    legend({'Target', 'Standard'},'Location','northwest')


 %pref-nonpref difference waves
erp_diff_out_face_in = squeeze(erp_out_face_in(:,1,:,:)-erp_out_face_in(:,2,:,:));
erp_diff_out_face_out = squeeze(erp_out_face_out(:,1,:,:)-erp_out_face_out(:,2,:,:));
 
subplot(1,3,3);
    boundedline(EEG.times,squeeze(mean(erp_diff_out_face_in(:,electrode,:),3)),squeeze(std(erp_diff_out_face_in(:,electrode,:),[],3))./sqrt(nsubs),'m',... %erp_diff_out_p(:,electrode,:),[],3
        EEG.times,squeeze(mean(erp_diff_out_face_out(:,electrode,:),3)),squeeze(std(erp_diff_out_face_in(:,electrode,:),[],3))./sqrt(nsubs),'b'); %erp_diff_out_np(:,2,electrode,i_cond,:),[],4   
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Difference-Wave by orientation');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    legend({'Facing Inside', 'Facing Outside'},'Location','northwest')
    
    
%Targets on same axis by facing orientation
figure;
boundedline(EEG.times,squeeze(mean(erp_out_face_in(:,1,electrode,:),4)),squeeze(std(erp_out_face_in(:,1,electrode,:),[],4))./sqrt(nsubs),'r',...
         EEG.times,squeeze(mean(erp_out_face_out(:,1,electrode,:),4)),squeeze(std(erp_out_face_out(:,1,electrode,:),[],4))./sqrt(nsubs),'b')

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
legend(L, {'Facing Inside', 'Facing Outside'},'Location','northeast')

    
%ttest for targets by preference at N1
time_window = find(EEG.times>125,1)-1:find(EEG.times>225,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_p(time_window,1,electrode,:),1)),squeeze(mean(erp_out_np(time_window,1,electrode,:),1)),.05,'both',1) 

%ttest for targets by preference at P2
time_window = find(EEG.times>225,1)-1:find(EEG.times>325,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_p(time_window,1,electrode,:),1)),squeeze(mean(erp_out_np(time_window,1,electrode,:),1)),.05,'both',1) 

%Standards on same axis
figure;
boundedline(EEG.times,squeeze(mean(erp_out_face_in(:,2,electrode,:),4)),squeeze(std(erp_out_face_in(:,2,electrode,:),[],4))./sqrt(nsubs),'r',...
         EEG.times,squeeze(mean(erp_out_face_out(:,2,electrode,:),4)),squeeze(std(erp_out_face_out(:,2,electrode,:),[],4))./sqrt(nsubs),'b')  
set(gca,'Color',[1 1 1]);
set(gca,'YDir','reverse');
axis tight; ylim([-8 12]);
line([-200 1000],[0 0],'color','k');
line([0 0],[-10 15],'color','k');
title('Standard Tones');
xlabel('Time (ms)');
ylabel('Voltage (uV)');

L(1) = plot(nan, nan, 'r');
L(2) = plot(nan, nan, 'b');
legend(L, {'Facing Inside', 'Facing Outside'},'Location','northeast')


%ttest for targets by preference at N1
time_window = find(EEG.times>125,1)-1:find(EEG.times>225,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_p(time_window,2,electrode,:),1)),squeeze(mean(erp_out_np(time_window,2,electrode,:),1)),.05,'both',1) 

%ttest for targets by preference at P2
time_window = find(EEG.times>225,1)-1:find(EEG.times>325,1)-2;
[h p ci stat] = ttest(squeeze(mean(erp_out_p(time_window,2,electrode,:),1)),squeeze(mean(erp_out_np(time_window,2,electrode,:),1)),.05,'both',1) 
    %%
%     %Averaging between 2 electrodes
%  electrode = [13, 14];
% figure
%     boundedline(EEG.times,squeeze(mean(mean(erp_out_face_in(:,1,electrode,:),3),4)),squeeze(std(erp_out_face_in(:,1,electrode,:),[],4))./sqrt(nsubs),'m',...
%         EEG.times,squeeze(mean(mean(erp_out_face_in(:,2,electrode,:),3),4)),squeeze(std(erp_out_face_in(:,2,electrode,:),[],4))./sqrt(nsubs),'k');    
%     set(gca,'Color',[1 1 1]);
%     set(gca,'YDir','reverse');
%     axis tight; ylim([-8 12]);
%     line([-200 1000],[0 0],'color','k');
%     line([0 0],[-2.5 8],'color','k');
%     title('Facing in');
%     xlabel('Time (ms)');
%     ylabel('Voltage (uV)');
%      legend({'Target', 'Standard'},'Location','northwest')
% 
%  
% figure
%     boundedline(EEG.times,squeeze(mean(mean(erp_out_face_out(:,1,electrode,:),3),4)),squeeze(std(erp_out_face_out(:,1,electrode,:),[],4))./sqrt(nsubs),'b',...
%         EEG.times,squeeze(mean(mean(erp_out_face_out(:,2,electrode,:),3),4)),squeeze(std(erp_out_face_out(:,2,electrode,i_cond,:),[],4))./sqrt(nsubs),'k');    
%     set(gca,'Color',[1 1 1]);
%     set(gca,'YDir','reverse');
%     axis tight; ylim([-8 12]);
%     line([-200 1000],[0 0],'color','k');
%     line([0 0],[-2.5 8],'color','k');
%     title('Facing out');
%     xlabel('Time (ms)');
%     ylabel('Voltage (uV)'); 
%     legend({'Target', 'Standard'},'Location','northwest')
% 
% 
%  %pref-nonpref difference waves
% erp_diff_out_face_in = squeeze(erp_out_face_in(:,1,:,:)-erp_out_face_in(:,2,:,:));
% erp_diff_out_face_out = squeeze(erp_out_face_out(:,1,:,:)-erp_out_face_out(:,2,:,:));
%  
% figure
% boundedline(EEG.times,squeeze(mean(mean(erp_diff_out_face_in(:,electrode,:),2),3)),squeeze(std(0))./sqrt(nsubs),'m',... %erp_diff_out_p(:,electrode,:),[],3
%         EEG.times,squeeze(mean (mean(erp_diff_out_face_out(:,electrode,:),2),3)),squeeze(std(0))./sqrt(nsubs),'b'); %erp_diff_out_np(:,2,electrode,i_cond,:),[],4   
%     set(gca,'Color',[1 1 1]);
%     set(gca,'YDir','reverse');
%     axis tight; ylim([-8 12]);
%     line([-200 1000],[0 0],'color','k');
%     line([0 0],[-2.5 8],'color','k');
%     title('Difference-Wave by orientation');
%     xlabel('Time (ms)');
%     ylabel('Voltage (uV)');
%     legend({'Facing Inside', 'Facing Outside'},'Location','northwest')
%     


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
boundedline(EEG.times,squeeze(mean(mean(erp_diff_out(:,electrode,:,:),3),4)), squeeze(std(0))./sqrt(nsubs),'m'),...
    
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


