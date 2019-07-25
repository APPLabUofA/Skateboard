%%
%%
ccc
exp = 'Skateboard';
subs = {'100' '102' '103' '104' '107' '108' '109' ...
    '111' '113' '116' '117' '118' '119' ...
    '125' '126'};
% subs = {'100' '101' '102' '103' '104' '106' '107' '108' '109' ...
%     '110' '111' '112' '113' '114' '115' '116' '117' '118' '119' ...
%     '120' '122' '123' '124' '125' '126' '127' '128' '129'};

nsubs = length(subs);
conds =  {'facing_In';'facing_Out'};
conds_lab = {'Facing Inside Track'; 'Facing Outside Track'};
nconds = length(conds);
Pathname = 'M:\Data\Skateboard\Winter2019\'; %M:\Data\Skateboard\Winter2019
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


for i_sub = 1:nsubs
    for i_cond = 1:nconds
        
        Filename = [subs{i_sub} '_' exp '_' conds{i_cond}];
        EEG = pop_loadset('filename',[Filename '_Corrected_Target.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_IO_V2\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        EEG = pop_loadset('filename',[Filename '_Corrected_Standard.set'],'filepath','M:\Data\Skateboard\Winter2019\segments_IO_V2\');
        [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
        
        
    end
end
eeglab redraw
%%
%subject erps
electrode = 15;%this is PZ in this electrode map
erp_out = [];
for i_sub = 1:nsubs
    fprintf(['Subject ' num2str(i_sub)])
    for i_cond = 1:nconds
        %average over trials (3rd dimension)
        erp_out(:,1,:,i_cond,i_sub) = mean(ALLEEG(1+ 2*((i_sub-1)*nconds+(i_cond-1))).data,3)'; %Targets
        erp_out(:,2,:,i_cond,i_sub) = mean(ALLEEG(2+ 2*((i_sub-1)*nconds+(i_cond-1))).data,3)'; %standards
    end
end

%%
%grand average plots + difference
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
figure
for i_cond = 1:nconds
    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'r';
    end
    
    subplot(2,nconds,i_cond);
    boundedline(EEG.times,squeeze(mean(erp_out(:,1,electrode,i_cond,:),5)),squeeze(std(erp_out(:,1,electrode,i_cond,:),[],5))./sqrt(nsubs),colour,...
        EEG.times,squeeze(mean(erp_out(:,2,electrode,i_cond,:),5)),squeeze(std(erp_out(:,2,electrode,i_cond,:),[],5))./sqrt(nsubs),'k');
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    if i_cond == 2
        legend('Targets','Standards','Location','NorthEast');
    end
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(conds_lab{i_cond});
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    
    subplot(2,nconds,nconds+i_cond);
    boundedline(EEG.times,squeeze(mean(erp_diff_out(:,electrode,i_cond,:),4)),squeeze(std(erp_diff_out(:,electrode,i_cond,:),[],4))./sqrt(nsubs),colour);
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    if i_cond == 2
        legend('Difference wave','Location','NorthEast');
    end
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title(conds_lab{i_cond});
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 DANIEL'S  PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Target data on same axis
figure
for i_cond = 1:nconds
    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'r';
    end
    
    boundedline(EEG.times,squeeze(mean(erp_out(:,1,electrode,i_cond,:),5)),squeeze(std(0))./sqrt(nsubs),colour); %std error erp_out(:,1,electrode,i_cond,:),[],5
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Targets by Orientation');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    
end
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r');
legend(L, {'facing in', 'facing out'})
hold on
%%
%%

%difference on same axis
erp_diff_out = squeeze(erp_out(:,1,:,:,:)-erp_out(:,2,:,:,:));
figure
for i_cond = 1:nconds
    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'r';
    end
    
    subplot;
    boundedline(EEG.times,squeeze(mean(erp_diff_out(:,electrode,i_cond,:),4)),squeeze(std(0))./sqrt(nsubs),colour);% std erp_diff_out(:,electrode,i_cond,:),[],4
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    axis tight; ylim([-8 12]);
    line([-200 1000],[0 0],'color','k');
    line([0 0],[-2.5 8],'color','k');
    title('Difference Wave By Facing Orientation');
    xlabel('Time (ms)');
    ylabel('Voltage (uV)');
    
end

L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r');
legend(L, {'dif. wave facing in', 'dif. wave facing out'})
hold on
%%
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


%%
%difference topographys
time_window = find(EEG.times>250,1)-1:find(EEG.times>450,1)-2;
figure('Color',[1 1 1]);
for i_cond = 1:nconds
    subplot(1,nconds,i_cond);
    set(gca,'Color',[1 1 1]);
    temp = mean(mean(erp_diff_out(time_window,:,i_cond,:),4),1)';
    temp(16:18) = NaN;
    topoplot(temp,'M:\Analysis\Skateboard\Skate_Vamp_Active_16.ced', 'whitebk','on','plotrad',.6,'maplimits',[-4 4]  )
    title(conds_lab{i_cond});
    t = colorbar('peer',gca);
    set(get(t,'ylabel'),'String', 'Voltage Difference (uV)');
    
end
%%

% for i_set = 1:48; trial_count(i_set) = ALLEEG(i_set).trials; end
% trial_count = reshape(trial_count,[2,3,8]);
% min(trial_count,[],3)
% mean(trial_count,3)
% max(trial_count,[],3)
%
% %mean and sd
% mean(mean(erp_diff_out(time_window,7,1:3,:),1),4)
% std(mean(erp_diff_out(time_window,7,1:3,:),1),[],4)
%
%
%
% ttest of each condition
[h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,7,1,:),1)),0,.05,'right',1)
[h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,7,2,:),1)),0,.05,'right',1)
[h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,7,3,:),1)),0,.05,'right',1)
[h p ci stat] = ttest(squeeze(mean(erp_diff_out(time_window,7,4,:),1)),0,.05,'right',1)



eeglab redraw
